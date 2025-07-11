args <- commandArgs(trailingOnly = TRUE)

#datafram1
bamstat2 <- read.table(args[1], header=T, comment.char = "")
#dataframe2
bamstat2 <- read.table(args[2], header=T, comment.char = "")
#metadata
samplemd <- read.table("../sample_data.txt", header=T)
