# import data
# Dataset mouse cortex GSE107122
library(data.table)


# Mouse 11.5
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2861nnn/GSM2861510/suppl/GSM2861510_E115_Only_Cortical_Cells_DGE.txt.gz",
              paste(data_dir,"E115_only_cortical_cells.txt.gz", sep = "/"),method = "wget", quiet = F)

# Mouse 13.5
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2861nnn/GSM2861511/suppl/GSM2861511_E135_Only_Cortical_Cells_DGE.txt.gz",
              paste(data_dir,"E135_only_cortical_cells.txt.gz", sep = "/"),method = "wget", quiet = F)


# Mouse 15.5
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107122/suppl/GSE107122_E155_Combined_Only_Cortical_Cells_DGE.txt.gz",
              paste(data_dir,"E155_only_cortical_cells.txt.gz", sep = "/"),method = "wget", quiet = F)

# Mouse 17.5
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2861nnn/GSM2861514/suppl/GSM2861514_E175_Only_Cortical_Cells_DGE.txt.gz",
              paste(data_dir,"E175_only_cortical_cells.txt.gz", sep = "/"),method = "wget", quiet = F)



raw = fread(paste(data_dir,"E175_only_cortical_cells.txt.gz", sep = "/"))
raw=as.data.frame(raw)
rownames(raw)=raw$V1
raw = raw[,2:ncol(raw)]

# Mouse dentate gyrus
# GSE104323
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104323/suppl/GSE104323_10X_expression_data_V2.tab.gz",
              paste(data_dir,"mouse_dentate_gyrus.tab.gz", sep = "/"),method = "wget", quiet = F)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104323/suppl/GSE104323_metadata_barcodes_24185cells.txt.gz",
              paste(data_dir,"metadata_mouse_dentate_gyrus.txt.gz", sep = "/"),method = "wget", quiet = F)



data = as.data.frame(fread(paste(data_dir,"mouse_dentate_gyrus.tab.gz", sep = "/"),sep = "\t"))
rownames(data) = data$cellid
metadata = read.csv(paste(data_dir,"metadata_mouse_dentate_gyrus.txt.gz", sep = "/"),header = T,sep = "\t")

#set the time point of interest and get, from the metadata, the cell's codes:
t = "E16.5"
cells.to.get = as.vector(metadata[metadata$characteristics..age == t,]$Sample.name..24185.single.cells.)


raw = data[, colnames(data) %in% cells.to.get]

if(!file.exists(paste(data_dir,"separation_age",sep = "/"))){   
  dir.create(file.path(paste(data_dir,"separation_age",sep = "/")))
}

write.csv(raw,paste(data_dir,"/separation_age/",t,"_mouse_dentate_gyrus.csv" ,sep = ""))



