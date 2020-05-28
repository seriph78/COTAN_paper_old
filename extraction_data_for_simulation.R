source("../src/functions.R")
t= "17"
data = load_files("../results/2019.12.05/",t = "17",s = "M1")# To load all data regarding the time point
raw =data[[1]]
cells=data[[2]] 
nu_est=data[[3]]
raw_norm = data[[4]] 
tot2 =data[[5]] 
mu_estimator = data[[6]] 
lambda_i = data[[7]]
dsq.e14ctx = data[[8]]


out_dir = "../results/2020.02.02/"
set =  unique(dsq.e14ctx@meta$clust)
#rownames(nu_est)= paste("E17_", rownames(nu_est), sep = "")
cluster_raw = nu_est
cluster_raw = cbind(cluster_raw,t(raw[rownames(cells),colnames(raw) %in% rownames(nu_est)]))
cluster_raw = t(cluster_raw)
rownames(cluster_raw)[1] = "nu"
for (c in set) {
  print(paste("cluster",c, sep=" "))
  n_set = c
  cells_set =  colnames(dsq.e14ctx@data[dsq.e14ctx@meta$clust == n_set] )
  if (length(cells_set) != length(unique(cells_set))) {
    print("error: duplicated cells!")
    break()
  }
  cluster_raw2 = cluster_raw[,colnames(cluster_raw) %in% cells_set]
  
  #write.csv(cluster_raw2, file=paste(out_dir,t,"_sqrt_cl_",c,".csv",sep = ""))
  
}
