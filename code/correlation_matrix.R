library(ggplot2)
library(reshape2)

data=read.table(file="DatasetM.txt", header=TRUE)
# Non-informative variable
data <- subset(data, select = -ID)
# Filter out censored data
com_data <- data[data$CEN == 1,] # Patients that have died

# Remove binnnary variables from dataset
cen_data <- data[, !(names(data) %in% c("HEP", "ASI", "DRUG", "CEN"))]
com_data <- com_data[, !(names(cen_data) %in% c("HEP", "ASI", "DRUG", "CEN"))]
#Get correlation matrix. Don't forget to ignore categorical variables in this process
cormat<-round(cor(cen_data),2)

################### AUXILIARY FUNCTIONS ##########################

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

##################################################################

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
{
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Pearson\nCorrelation") +
    theme_minimal() + # minimal theme
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, 
                                 size = 12, hjust = 1),
      legend.text = element_text(size = 12),  # Increase legend text size
      legend.title = element_text(size = 12)  # Increase legend title size
    ) +
    coord_fixed()
  
  # Adds correlation coefficients on the heatmap
  ggheatmap <- ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 6) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
  print(ggheatmap)
}


ggsave(file="corr_matrix.pdf", width=40, height=40, dpi=300)


#Remove one variable of pair with high correlation.
#Do correlation with response variable (watch out because response variable is categorical)











# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Create a ggheatmap
{
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", guide = "none") +
    ggtitle("Correlation Heatmap\nfor all patients") +  # Add a title to the plot
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, 
                                 size = 12, hjust = 1),
      legend.text = element_text(size = 12),  # Increase legend title size
      plot.title = element_text(size = 16)  # Increase title size
    ) +
    coord_fixed()
  
  # Adds correlation coefficients on the heatmap
  ggheatmap <- ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 6) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"  # Remove the legend
    )
  
  print(ggheatmap)
}



