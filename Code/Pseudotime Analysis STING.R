library(zellkonverter)

url <- "https://datasets.cellxgene.cziscience.com/0d9f6c76-50d2-42ab-bcb7-10ab2847d2fd.h5ad"
tmp_file <- tempfile(fileext = ".h5ad")
download.file(url, destfile = tmp_file, mode = "wb")
