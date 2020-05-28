library(parallel)
library(data.table)

nu.samples = read.csv("../results/2020.01.31/17_sqrt/nu_out17_sqrt.csv", header = T,row.names = 1)
cl = unique(nu.samples$cluster)
cl.perc = vector(length = length(cl))
names(cl.perc) = cl
for (c in cl) {
   cl.perc[c] = length(nu.samples$cluster[nu.samples$cluster == c])/length(nu.samples$cluster)
}
sum(cl.perc)
 
eta = read.csv("../results/2020.01.31/17_sqrt/eta_out17_sqrt.csv", header = T,row.names = 1,na.strings=c("","NA"), stringsAsFactors = F)
colnames(eta) = eta[1,]
eta = eta[2:nrow(eta),]
eta=data.matrix(eta)



theta = read.csv("../results/2020.01.31/17_sqrt/theta_out17_sqrt.csv", header = T,row.names = 1,na.strings=c("","NA"), stringsAsFactors = F)
colnames(theta) = theta[1,]
theta = theta[2:nrow(theta),]
theta=data.matrix(theta)

#-----------------------------------
# not updated
generate.cell.single = function(c){
  nu.val = sample(nu.samples$nu,1)
  
  #cell.group = 1
  cell.group = sample(c(1,2),1) #to create a population of amizure of the first two clusters
  
  cells = data.frame(rep(NA, (nrow(eta))),row.names = rownames(eta)[1:nrow(eta)])
  rownames(cells)[1]="nu"
  for(g in rownames(cells)){
    if (g == "nu") {
      cells[g,1] =nu.val
    }else{
      cells[g,1] =  rpois(1,rgamma(1,shape = eta[g,cell.group], scale = nu.val*theta[g,cell.group]))
    }
  }
  colnames(cells) = paste("cell",c,sep = "_")
  return(cells)
}

#-------------------------------------

generate.cell = function(c){
  nu.val = sample(nu.samples$nu,1)
  
  r = runif(1,0,1)
  #cell.group = names(cl.perc)[sum(r >= cumsum(cl.perc))+1]
  cell.group = names(which(cl.perc == max(cl.perc))) # to generate a single cluster dataset (using the biggest)
  
  cells = data.frame(rep(NA, (nrow(eta))),row.names = rownames(eta))
  rownames(cells)[1]="nu"
  for(g in rownames(cells)){
    if (g == "nu") {
      cells[g,1] =nu.val
    }else{
      cells[g,1] =  rpois(1,rgamma(1,shape = eta[g,cell.group], scale = nu.val*theta[g,cell.group]))
      }
    }
  colnames(cells) = paste("cell",c,cell.group,sep = "_")
  return(cells)
}

cells.sim.base = mclapply(c(1:4000), generate.cell,mc.cores = 11)

cells.sim = cells.sim.base[[1]]
for (tt in 2:length(cells.sim.base)) {
  cells.sim= cbind(cells.sim,cells.sim.base[[tt]])
  
}

rm(cells.sim.base)
gc()

fwrite(cells.sim, "../results/2020.02.02/simulated.dataset_E17_sqrt.fra.1cl_4000cs.csv.gz",compress = "gzip", row.names = T, col.names = T)


