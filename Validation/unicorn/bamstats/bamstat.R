args <- commandArgs(trailingOnly = TRUE)

#datafram1
bamstat1 <- read.table(args[1], header=T, comment.char = "")
#dataframe2
bamstat2 <- read.table(args[2], header=T, comment.char = "")
#metadata
samplemd <- read.table("../sample_data.txt", header=T)

#merge data by existing samples
merged <- merge(bamstat1, bamstat2, by = "X.name", suffixes = c(".1", ".2"))
# Add site and library type info by matching X.name to library_id in samplemd
merged$site <- samplemd$site[match(merged$X.name, samplemd$library_id)]
merged$library_type <- samplemd$library_type[match(merged$X.name, samplemd$library_id)]
site_levels  <- unique(merged$site)
site_colors  <- setNames(rainbow(length(site_levels)), site_levels)
point_colors <- site_colors[merged$site]
ltype_levels <- unique(merged$library_type)
ltype_shapes <- setNames(c(19,17), ltype_levels)
point_shapes <- ltype_shapes[merged$library_type]

#plot
svg(args[3], width = 8, height = 6) # Save plot to SVG file
plot(
  merged$mani.1, merged$mani.2,
  xlab = "Ref 1.0",
  ylab = "Ref 3.0",
  main = "mean ANI: Ref1.0 vs Ref3.0",
  pch = point_shapes, col = point_colors
)
abline(0, 1, col = "red", lty = 2) # Optional: reference line y=x

leg1<-legend("topleft", legend = site_levels,
       col = site_colors, pch = 19, title = "Site")
legend(leg1$rect$left + leg1$rect$w, 
       leg1$rect$top, legend = ltype_levels,
       pch = ltype_shapes, title = "Library Type")
dev.off() # Close the SVG device