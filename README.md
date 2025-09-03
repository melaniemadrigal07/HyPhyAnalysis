# HyPhy Analysis Toolkit
Post-HyPhy pipeline for extracting and analyzing fungal hyphal network data, which can be used after running: https://github.com/AaronMoseley/SkeletonizationTool.git.

This repository houses two core scripts:
1. **countspores.R** - Counts spores per image using a fixed grid, then estimates **spores per nL** and **spores per µL**.
2. **HyPhyPipeline.R** - extracts metadata and network traits from skeletonized images **and plots networks** by connecting points and building edges.


## Requirements
- R (≥ 4.2.0)  
- Packages: `tidyverse`, `igraph`, `vroom`, `ggplot2`

## Inputs
1. **countspores.R**  
   - A folder containing **PNG images** of spores.  
   - The script processes each image, applies filtering and thresholding, then counts spores per grid square.
   - Outputs a .csv file with spore counts of each sample, which can be used for further dilution. 

2. **HyPhyPipeline.R**  
   - Each subfolder ending in `skeleton_csvs/` must include:  
     - `fileInfo.csv`  
     - `network_points.csv`  
     - `network_lines.csv`  
     - `network_metadata.csv`
       
## Installation
To use this repository, clone it by:
```bash
git clone https://github.com/yourname/HyPhyAnalysis.git
cd HyPhyAnalysis
