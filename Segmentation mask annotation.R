# This code performs single cell annotation on images with follicular B cell segmentation mask created by Cellpose

install.packages("BiocManager") 
BiocManager::install("EBImage")

library(ggplot2)
library(readr)
library(dplyr)
library(imager)
library(png)
library(tiff)
library(EBImage)

# load region properties tables for CD45.1 and CD45.2 channels
CD45_1_props <- read_csv("Cellpose_training/Napari_output/CTransfer time course R1_D5pt_spl1_3 30x_1 stitch merge-MZp7 -CD169_CD45.1RegionPropsTable.csv")
CD45_2_props <- read_csv("Cellpose_training/Napari_output/CTransfer time course R1_D5pt_spl1_3 30x_1 stitch merge-MZp7 -CD169_CD45.2RegionPropsTable.csv")
B220_props <- read_csv("Cellpose_training/Napari_output/CTransfer time course R1_D5pt_spl1_3 30x_1 stitch merge-MZp7 -CD169_B220RegionPropsTable.csv")

# merge the tables based on the "label" column
merged_props <- CD45_1_props %>% 
  inner_join(CD45_2_props, by = "label", suffix = c("_CD45_1", "_CD45_2"))

# plot histograms of mean intensities for both channels 
ggplot(merged_props, aes(x = mean_intensity_CD45_1)) + 
  geom_histogram(binwidth = 10) + 
  labs(title = "Distribution of Mean Intensity CD45.1", x = "Mean Intensity CD45.1", y = "Frequency")

ggplot(merged_props, aes(x = mean_intensity_CD45_2)) + 
  geom_histogram(binwidth = 10) + 
  labs(title = "Distribution of Mean Intensity CD45.2", x = "Mean Intensity CD45.2", y = "Frequency")

# make scatter plot to visualize relationship between the intensities of the two channels
ggplot(merged_props, aes(x = mean_intensity_CD45_1, y = mean_intensity_CD45_2)) +
  geom_point() +
  labs(title = "Mean Intensity Scatter Plot", x = "Mean Intensity CD45.1", y = "Mean Intensity  CD45.2")

# Determine the range of the intensity column
range_CD45_1 <- range(merged_props$mean_intensity_CD45_1, na.rm = TRUE)
range_CD45_2 <- range(merged_props$mean_intensity_CD45_2, na.rm = TRUE)

# Generate histogram data
hist_CD45_1 <- hist(merged_props$mean_intensity_CD45_1, breaks = seq(floor(range_CD45_1[1]), ceiling(range_CD45_1[2]), by = 1), plot = TRUE)
hist_CD45_2 <- hist(merged_props$mean_intensity_CD45_2, breaks = seq(floor(range_CD45_2[1]), ceiling(range_CD45_2[2]), by = 1), plot = TRUE)

# Create an image-like structure (2D matrix) from the histogram counts
hist_image_CD45_1 <- matrix(hist_CD45_1$counts, nrow = length(hist_CD45_1$counts))
hist_image_CD45_2 <- matrix(hist_CD45_2$counts, nrow = length(hist_CD45_2$counts))

# Set intensity thresholds for positive identification
threshold_CD45_1_value <- 1000
threshold_CD45_2_value <- 1000

# Define your annotation logic with the intensity thresholds
cell_annotation <- function (mean_intensity_CD45_1, mean_intensity_CD45_2) {
  if (mean_intensity_CD45_1 >= threshold_CD45_1_value & mean_intensity_CD45_2 >= threshold_CD45_2_value)
  {
    return("WT donor")
  } else if (mean_intensity_CD45_1 >= threshold_CD45_1_value & mean_intensity_CD45_2 < threshold_CD45_2_value) {
    return ("564 recipient") 
  } else if (mean_intensity_CD45_1 < threshold_CD45_1_value & mean_intensity_CD45_2 >= threshold_CD45_2_value) {
    return ("TLR7KO donor")
  } else {
    return("unclassified")
  }
}

# Apply the annotation logic
merged_props <- merged_props %>% 
  rowwise() %>% 
  mutate(annotation = cell_annotation(mean_intensity_CD45_1, mean_intensity_CD45_2))

# Save the annotated table
# write_csv(merged_props, "MZ_annotated_cells.csv")

######### Create annotated label image ########

# Load the segmentation mask
segmentation_mask <- load.image("Cellpose_training/CTransfer_time_course_-CD169_segmentation/CTransfer time course R1_D5pt_spl1_3 30x_1 stitch merge-MZp7 -CD169_cp_masks.png")

# Convert the segmentation mask to a matrix
segmentation_matrix <- as.matrix(segmentation_mask)

# Scale up values in matrix based on 16 bit image largest values
segmentation_matrix <- segmentation_matrix * 65535

# Create an empty matrix for the annotated labels 
annotated_label_matrix <- matrix(0, nrow = nrow(segmentation_matrix), ncol = ncol(segmentation_matrix))

# Map the annotations to the label matrix 
annotation_mapping <- list("WT donor" = 10, "TLR7KO donor" = 20, "564 recipient" = 30, "unclassified" = 0)

for(i in 1:nrow(merged_props)) {
  label <- merged_props$label[i]
  annotation <- merged_props$annotation[i] 
  
  # Check if any elements in segmentation_matrix match the current label 
  matching_indices <- which(segmentation_matrix == label & segmentation_matrix != 0)
  
  print(paste("Processing label:", label, "with annotation:", annotation))
  print(paste("Number of matching elements:", length(matching_indices)))
  print(paste("Matching indices:", matching_indices))
  
  if (length(matching_indices) > 0) {
    # Perform the assignment only if there are matching elements
    annotated_label_matrix[segmentation_matrix == label] <- annotation_mapping[[annotation]]
  } else {
    print(paste("No elements in segmentation_matrix match label:", label))
  }
}

# Save the annotated label image 
writePNG(annotated_label_matrix, "Cellpose_training/CTransfer_time_course_-CD169_segmentation/Annotated masks/CTransfer time course R1_D5pt_spl1_3 30x_1 stitch merge-MZp7 -CD169_annotated_segmentation_mask.png")






