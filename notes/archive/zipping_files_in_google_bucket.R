# zipping files in Google Bucket

library(googleCloudStorageR)

gcs_auth()
gcs_global_bucket("s05-ad-sample-data")

files <- gcs_list_objects(prefix = "RNAseq_24hr/")[-1,]


for (file in files$name) {
  # Download file
  gcs_get_object(file, saveToDisk = basename(file), overwrite = TRUE)
  
  # Gzip the file
  gz_file <- paste0(basename(file), ".gz")
  R.utils::gzip(basename(file), destname = gz_file, overwrite = TRUE)
  
  # Upload gzipped file
  gcs_upload(gz_file, name = paste0(file, ".gz"))
  
  # Clean up
  file.remove(basename(file), gz_file)
  gcs_delete_object(file)
}
