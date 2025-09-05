if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("EBImage")
library(EBImage)

# Parameters (make sure to change image_dir to where your images are located)
image_dir <- "/Users/melaniemadrigal/Desktop/dowell/spore\ counting/peter/peter/8.11.25"  # Change to your directory
image_files <- list.files(image_dir, pattern = "\\.png$", full.names = TRUE)

#count paramaters optimized for Cytation 5
n_rows <- 5  # Number of rows in the grid

n_cols <- 5  # Number of columns in the grid
volume_per_square_nl <- 50 
pixels_per_mm <- 2850  # Pixels per millimeter (for scaling purposes) (obtained from imagej)
mm_per_square <- 0.2  # Side length of each square in mm (for grid)

# ---- Results Container ----
results <- data.frame(
  Image = character(),
  AverageSpores = numeric(),
  Concentration_Sp_nL = numeric(),
  Concentration_Sp_uL = numeric(),
  stringsAsFactors = FALSE
)

# ---- Loop Over Images ----
for (image_path in image_files) {
  cat("Processing:", basename(image_path), "\n")
  
  # Load and preprocess the image
  img <- readImage(image_path)
  gray <- channel(img, "gray")  # Convert to grayscale
  norm <- normalize(gray)  # Normalize the image
  blurred <- gblur(norm, sigma = 1)  # Apply Gaussian blur to smooth the image
  
  # Thresholding to create a binary image with adjusted parameters
  binary <- thresh(blurred, w=10, h=10, offset=0.02)  # Adjusting the window size and offset
  
  # Apply closing to join together objects that are part of the same spore, if needed
  closed <- closing(binary, makeBrush(5, shape='disc'))
  
  # Apply opening to remove small noise objects
  cleaned <- opening(closed, makeBrush(5, shape='disc'))
  
  # Label connected components (objects/spores)
  labeled <- bwlabel(cleaned)
  
  # ---- Get Image Dimensions and Grid Size ----
  dims <- dim(labeled)
  box_h <- dims[1] / n_rows  # Height of each grid box
  box_w <- dims[2] / n_cols  # Width of each grid box
  
  # ---- Count Spores in Each Grid Square ----
  box_counts <- numeric()  # Container for spore counts in each grid box
  box_id <- 1  # Grid box identifier
  
  # Count spores in each grid square
  for (i in 0:(n_rows - 1)) {
    for (j in 0:(n_cols - 1)) {
      x1 <- as.integer(j * box_w) + 1
      x2 <- as.integer((j + 1) * box_w)
      y1 <- as.integer(i * box_h) + 1
      y2 <- as.integer((i + 1) * box_h)
      
      # Skip grid boxes that go out of image bounds
      if (x2 > dims[2] || y2 > dims[1]) next
      
      # Extract the sub-image corresponding to the grid box
      sub_img <- labeled[y1:y2, x1:x2]
      
      # Count the number of connected components (spores) in the sub-image
      spore_count <- max(bwlabel(sub_img))  # max of bwlabel gives the number of objects
      box_counts[box_id] <- spore_count  # Store the count for this grid box
      box_id <- box_id + 1
    }
  }
  
  # Calculate Metrics
  average_spores <- mean(box_counts, na.rm = TRUE)  # Average spore count across grid boxes
  concentration_nl <- average_spores / volume_per_square_nl  # Concentration in spores per nanoliter
  concentration_ul <- concentration_nl * 1000  # Convert to spores per microliter
  
  # Store results in the results dataframe
  results <- rbind(results, data.frame(
    Image = basename(image_path),
    AverageSpores = average_spores,
    Concentration_Sp_nL = concentration_nl,
    Concentration_Sp_uL = concentration_ul,
    stringsAsFactors = FALSE
  ))
  
  # Display the processed image (after processing each image,visual check)
  display(colorLabels(labeled), method = "raster")
}

#  Output Results
print(results)
results$Concentration_Sp_uL
write.csv(df$Concentration_Sp_uL, "Concentration_Sp_uL.csv", 
            col.names = FALSE, row.names = FALSE, sep = ",")
#these results are then inputted into an excel sheet for diluations



