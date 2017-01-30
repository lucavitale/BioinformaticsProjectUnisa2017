files <- as.data.frame(list.files("BreastMaleData"))
fm <- read.table("file_matching.txt", header = TRUE, row.names = NULL, sep = "\t")
names(files) <- "file_id"
merged <- merge(files, fm, by = "file_id")
