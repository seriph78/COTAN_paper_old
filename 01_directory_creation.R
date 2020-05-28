# subdirectory formation
main_dir = "../rstudio/storage/cotan"
data_dir = paste(main_dir,"data",sep = "/")
res_dir = paste(main_dir,"results",sep = "/")

if(!file.exists(main_dir)){
  dir.create(file.path(main_dir))
}

if(!file.exists(data_dir)){   
  dir.create(file.path(data_dir))
}

if(!file.exists(res_dir)){   
  dir.create(file.path(res_dir))
}

setwd("/home/rstudio/")
