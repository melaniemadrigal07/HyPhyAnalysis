#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggplot2)
  library(stringr)
  library(wesanderson)
})

# EDIT this to where you have your data stored
base_dir <- "/Volumes/Lexar/HyphalNetDS/All"

# Map well IDs to isolate names (edit as needed)
well_map <- c(
  "C5"  = "B05.10",
  "G12" = "KatieTomato",
  "H6"  = "1.04.25"
)

min_tp_keep <- 4  # keep timepoints >= this

#helper functions
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

num_or_na <- function(x) suppressWarnings(readr::parse_number(as.character(x)))

# looking for folders with information to be extracted
# Find all folders that end with 'skeleton_csvs'
skeleton_folders <- Filter(
  function(x) grepl("skeleton_csvs$", x),
  list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
)

# Output dirs
dir.create("plots", showWarnings = FALSE)
dir.create("stats", showWarnings = FALSE)

# Storage
overall_stats <- tibble()
meta_traits   <- tibble()

# Main loop responsible for pulling info from paths
for (folder in skeleton_folders) {
  cat("Processing:", folder, "\n")
  
  # Paths
  fileinfo_path <- file.path(folder, "fileInfo.csv")
  points_path   <- file.path(folder, "network_points.csv")
  lines_path    <- file.path(folder, "network_lines.csv")
  md_path       <- file.path(folder, "network_metadata.csv")
  
  # Guards
  if (!file.exists(fileinfo_path)) { warning("Missing fileInfo.csv in: ", folder); next }
  if (!file.exists(points_path) || !file.exists(lines_path)) {
    warning("Missing points/lines CSV in: ", folder); next
  }
  
  # --- sample & tp from fileInfo.csv (read once) ---
  fi <- readr::read_csv(fileinfo_path, col_names = FALSE, show_col_types = FALSE)
  if (ncol(fi) < 2) { warning("Bad fileInfo.csv in: ", folder); next }
  colnames(fi) <- c("name", "value")
  
  sample_val <- fi %>% dplyr::filter(name == "sample")    %>% dplyr::pull(value)
  tp_raw     <- fi %>% dplyr::filter(name == "timestamp") %>% dplyr::pull(value)
  tp <- ifelse(is.na(tp_raw) || tp_raw == "",
               stringr::str_extract(fi$value[1], "\\d+"),
               tp_raw) %>% as.integer()
  
  # --- metadata traits (optional per folder) ---
  if (file.exists(md_path)) {
    md <- readr::read_csv(md_path, show_col_types = FALSE)
    
    md_keep <- md %>%
      dplyr::filter(name %in% c("fractalDimension", "linesInImage", "clustersInImage")) %>%
      dplyr::mutate(value_num = num_or_na(value)) %>%
      dplyr::select(name, value_num) %>%
      tidyr::pivot_wider(names_from = name, values_from = value_num) %>%
      dplyr::mutate(sample = sample_val, tp = tp)
    
    meta_traits <- dplyr::bind_rows(meta_traits, md_keep)
  } else {
    warning("Missing network_metadata.csv in: ", folder)
  }
  
  #  Load & expand lines (reconstructing hyphal segments)
  network_lines <- readr::read_csv(lines_path, col_names = FALSE, show_col_types = FALSE)
  if (ncol(network_lines) < 2) { warning("Bad network_lines.csv in: ", folder); next }
  colnames(network_lines)[1:2] <- c("line_id", "pointIndices")
  
  lines_long <- network_lines %>%
    mutate(line_id = suppressWarnings(as.integer(line_id)),
           pointIndices = str_replace_all(pointIndices, "[^0-9,]", "")) %>%
    filter(!is.na(line_id), !is.na(pointIndices), nchar(pointIndices) > 0) %>%
    separate_rows(pointIndices, sep = ",") %>%
    mutate(pointIndex = suppressWarnings(as.integer(str_trim(pointIndices)))) %>%
    filter(!is.na(pointIndex)) %>%
    group_by(line_id) %>%
    mutate(order = dplyr::row_number()) %>%
    ungroup()
  
  # Load points
  points <- readr::read_csv(points_path, show_col_types = FALSE) %>%
    mutate(pointIndex = as.integer(pointIndex))
  
  #Join path coords
  line_paths <- lines_long %>%
    left_join(points, by = "pointIndex") %>%
    filter(!is.na(x), !is.na(y)) %>%
    mutate(line_id = as.integer(line_id))
  
  if (nrow(line_paths) == 0) { warning("No line paths in: ", folder); next }
  
  # Build edges for graph
  global_edges <- line_paths %>%
    arrange(line_id, order) %>%
    group_by(line_id) %>%
    mutate(next_point = dplyr::lead(pointIndex)) %>%
    filter(!is.na(next_point)) %>%
    ungroup() %>%
    select(from = pointIndex, to = next_point)
  
  if (nrow(global_edges) == 0) { warning("No edges in: ", folder); next }
  
  g <- igraph::graph_from_data_frame(global_edges, directed = FALSE)
  
  # Node degree types
  degree_df <- tibble(
    pointIndex = as.integer(names(igraph::degree(g))),
    degree = as.integer(igraph::degree(g))
  )
  
  points_deg <- points %>%
    left_join(degree_df, by = "pointIndex") %>%
    mutate(degree = dplyr::coalesce(degree, 0L),
           type = case_when(
             degree == 1 ~ "Tip",
             degree >  2 ~ "Branch",
             degree == 0 ~ "Isolated",
             TRUE        ~ "Internal"
           ))
  
  # Stats from graph
  num_edges    <- igraph::gsize(g)
  num_tips     <- sum(points_deg$type == "Tip", na.rm = TRUE)
  num_branches <- sum(points_deg$type == "Branch", na.rm = TRUE)
  num_isolated <- sum(points_deg$type == "Isolated", na.rm = TRUE)
  
  cat("  Edges:", num_edges, " Tips:", num_tips,
      " Branches:", num_branches, " Isolated:", num_isolated, "\n")
  
  # --- Save per-folder stats row ---
  stats_row <- tibble(
    folder   = basename(folder),
    sample   = sample_val,
    timepoint = tp,
    edges    = num_edges,
    tips     = num_tips,
    branches = num_branches,
    isolated = num_isolated
  )
  overall_stats <- bind_rows(overall_stats, stats_row)
  
  # Plot & save network image
  p <- ggplot() +
    geom_path(data = line_paths %>% arrange(line_id, order),
              aes(x = x, y = y, group = factor(line_id)),
              color = "white", linewidth = 0.2, alpha = 0.6) +
    geom_point(data = points_deg, aes(x = x, y = y, color = type), size = 0.5) +
    scale_color_manual(values = c(
      "Tip" = "yellow",
      "Branch" = "limegreen",
      "Internal" = "deeppink",
      "Isolated" = "dodgerblue"
    )) +
    coord_equal() +
    theme_void(base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(color = "white"),
      plot.title = element_text(color = "white", size = 14, face = "bold")
    ) +
    labs(title = paste0("Hyphal Skeleton\nSample: ", sample_val, " | Timepoint: ", tp))
  
  ggsave(
    filename = paste0(basename(folder), "_network.png"),
    plot = p,
    path = "plots",
    width = 6, height = 6, dpi = 300
  )
}

# Write raw tables
readr::write_csv(overall_stats, "stats/network_summary_stats.csv")
readr::write_csv(meta_traits,  "stats/network_meta_traits_raw.csv")

#Summarize traits & join 
# 1) Metadata traits → isolate × tp (mean ± SE)
meta_summary <- meta_traits %>%
  mutate(
    isolate = dplyr::recode(sample, !!!well_map),
    tp = as.integer(tp)
  ) %>%
  filter(!is.na(isolate), !is.na(tp), tp >= min_tp_keep) %>%
  group_by(isolate, tp) %>%
  summarise(
    mean_fd = mean(fractalDimension, na.rm = TRUE),
    se_fd   = sd(fractalDimension,   na.rm = TRUE) / sqrt(sum(!is.na(fractalDimension))),
    mean_lI = mean(linesInImage,     na.rm = TRUE),
    se_lI   = sd(linesInImage,       na.rm = TRUE) / sqrt(sum(!is.na(linesInImage))),
    mean_cI = mean(clustersInImage,  na.rm = TRUE),
    se_cI   = sd(clustersInImage,    na.rm = TRUE) / sqrt(sum(!is.na(clustersInImage))),
    .groups = "drop"
  )

# 2) Hyphal stats → isolate × tp (mean ± SE)
hypha_summary <- overall_stats %>%
  mutate(
    isolate = dplyr::recode(sample, !!!well_map),
    tp = as.integer(timepoint)
  ) %>%
  filter(!is.na(isolate), !is.na(tp), tp >= min_tp_keep) %>%
  group_by(isolate, tp) %>%
  summarise(
    mean_tips     = mean(tips,     na.rm = TRUE),
    se_tips       = sd(tips,       na.rm = TRUE) / sqrt(sum(!is.na(tips))),
    mean_edges    = mean(edges,    na.rm = TRUE),
    se_edges      = sd(edges,      na.rm = TRUE) / sqrt(sum(!is.na(edges))),
    mean_branches = mean(branches, na.rm = TRUE),
    se_branches   = sd(branches,   na.rm = TRUE) / sqrt(sum(!is.na(branches))),
    .groups = "drop"
  )

# 3) Join
combo <- meta_summary %>% left_join(hypha_summary, by = c("isolate","tp"))

readr::write_csv(combo, "stats/isolate_tp_summary_fd_lines_clusters_tips_edges_branches.csv")

# quick quick to ensure plots make sense
# A) FD vs Tips
pal <- wes_palette("FantasticFox1", 5)
ggplot(combo, aes(mean_fd, mean_tips, color = isolate)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("B05.10"=pal[1], "KatieTomato"=pal[3], "1.04.25"=pal[2])) +
  labs(x = "Fractal dimension (mean)", y = "Tips (mean)", color = "Isolate",
       title = "FD vs Tips across timepoints") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.line = element_line()) +
  ggsave("plots/fd_vs_tips.png", width = 7, height = 5, dpi = 300)

# B) Standardized time-course per isolate
combo_long <- combo %>%
  select(isolate, tp, mean_fd, mean_lI, mean_cI, mean_tips, mean_edges, mean_branches) %>%
  pivot_longer(-c(isolate, tp), names_to = "metric", values_to = "value") %>%
  group_by(isolate, metric) %>%
  mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
  ungroup()

metric_cols <- c(
  mean_fd        = pal[5],  # red
  mean_lI        = pal[2],  # blue
  mean_cI        = pal[3],  # yellow
  mean_tips      = pal[1],  # orange
  mean_edges     = pal[4],  # gray
  mean_branches  = "#6A3D9A" # purple-ish
)

ggplot(combo_long, aes(tp, z, color = metric, group = metric)) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2) +
  facet_wrap(~ isolate, nrow = 1, scales = "free_x") +
  scale_color_manual(values = metric_cols, name = NULL) +
  scale_x_continuous(breaks = sort(unique(combo_long$tp))) +
  labs(x = "Timepoint", y = "Standardized (z-score)",
       title = "FD, Lines, Clusters & Hyphal traits over time (per isolate)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.line = element_line(),
        legend.position = "bottom") +
  ggsave("plots/standardized_timecourse.png", width = 10, height = 4.5, dpi = 300)

cat("Done. Wrote:\n",
    "- stats/network_summary_stats.csv (per-folder hyphal stats)\n",
    "- stats/network_meta_traits_raw.csv (per-folder metadata traits)\n",
    "- stats/isolate_tp_summary_fd_lines_clusters_tips_edges_branches.csv (joined summaries)\n",
    "Plus plots in ./plots\n")
