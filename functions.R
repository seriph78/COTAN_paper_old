library(dplyr)
library(reticulate)
library(propagate)
library(parallel)
library(data.table)
library(matrixStats)
library(rray)
library(tibble)
library(Seurat)

library(ggfortify)
library(ggplot2)
library(ggsci)
library(gmodels)
#-------------------------------------------

c1 = 6.5
c2 = 60   
## psi functions
psi <- function(x){
  val_psi <- x**2 + psi_1(x) * psi_2(x) + psi_3(x)
  return(val_psi)
}

psi_1 <- function(x){
  val_psi_1 = 1/4*x^7 + 3/32*x^5 + x - x^2/sqrt(2)
  return(val_psi_1)
}

psi_2 <- function(x){
  val_psi_2 = 1 / (x^7 + 1)
  return(val_psi_2)
}

psi_3 <- function(x){
  val_psi_3 = c2 * x^8 * exp(-c1*x)
  return(val_psi_3)
}

## derivates of psi
first_der_psi <- function(x){
  f_psi = 2 * x + first_der_psi_1(x) * psi_2(x) + psi_1(x) * first_der_psi_2(x) + first_der_psi_3(x)
  return(f_psi)
}

sec_der_psi <- function(x){
  s_psi = 2 + sec_der_psi_1(x) * psi_2(x) + 2 * first_der_psi_1(x) * first_der_psi_2(x) + psi_1(x) * sec_der_psi_2(x) + sec_der_psi_3(x)
  return(s_psi)
}

first_der_psi_1 <- function(x){
  f_psi_1 = 7/4 * x^6 + 15/32 * x^4 + 1 - sqrt(2) * x
  return(f_psi_1)
}

sec_der_psi_1 <- function(x){
  s_psi_1 = 21/2 * x^5 + 15/8 * x^3 - sqrt(2)
  return(s_psi_1)
}

first_der_psi_2 <- function(x){
  f_psi_2 = -(7 * x^6)/(x^7 + 1)^2
  return(f_psi_2)
}

sec_der_psi_2 <- function(x){
  s_psi_2 = -(14 * x^5 * (3 - 4*x^7))/(x^7 + 1)^3
  return(s_psi_2)
}

first_der_psi_3 <- function(x){
  f_psi_3 = c2*(8*x^7 - c1*x^8)*exp(-c1*x)
  return(f_psi_3)
}

sec_der_psi_3 <- function(x){
  s_psi_3 = c2 * (56 * x^6 - 16 * c1 * x^7 + c1^2 * x^8) * exp(-c1*x)
  return(s_psi_3)
}

# tau
tau <- function(x){
  xtau = psi(x) - x^2
  return(xtau)
}

# lambda
lambda <-function(xi_star, mat_X){
  v = rowVarsC(mat_X)
  l_i = psi(xi_star) + 1/2 * sec_der_psi(xi_star)*(v - tau(xi_star))
  return(l_i)
}

# nu
R_star_j <- function(x_star_j, mat_X){
  v = colVarsC(mat_X)
  Rsj = psi(x_star_j) + 1/2 * sec_der_psi(x_star_j)*(v - tau(x_star_j))
  return(Rsj)
}

nu <- function(x_star_j,mat_X){
  m = length(colnames(mat_X))
  Rsj = R_star_j(x_star_j, mat_X)
  nu_j = Rsj/(sum(Rsj)/m)
  return(nu_j)
}

#---------------------------------------------
fun_pzero <- function(a,mu){
  #a= as.numeric(a[,2])
  #print(a)
  (a <= 0)*(exp(-(1+abs(a))*mu)) + (a > 0)*(1+abs(a)*mu)^(-1/abs(a))
}

fun_pzero_posi <- function(r,mu){ (1+r*mu)^(-1/r) }

fun_pzero_nega0 <- function(r,mu){ (exp(-(1-r)*mu))}

fun_dif_mu_zeros <- function(h,x,somma_zeri){
  if (h > 0) {
    sum(fun_pzero_posi(h,mu_estimator[x,])) - somma_zeri#/somma_zeri  
  }else{
    sum(fun_pzero_nega0(h,mu_estimator[x,])) - somma_zeri#/somma_zeri  
  }
}

fun_my_opt <- function(x){
  somma_zeri = sum(cells[x,] == 0)
  #somma_zeri = rowSums(cells[x,] == 0)
  a1 = 0
  u1 = fun_dif_mu_zeros(a1,x,somma_zeri)
  a2 = a1
  u2 = u1
  if (u1 > 0) {
    a1 = a1 - 1
    u1 = fun_dif_mu_zeros(a1,x,somma_zeri)
    while (u1 > 0) {
      a2 = a1
      u2 = u1
      a1 = 2 * a1
      u1 = fun_dif_mu_zeros(a1,x,somma_zeri)
    }
  }else{
    a2 = 1
    u2 = fun_dif_mu_zeros(a2,x,somma_zeri)
    while (u2 < 0) {
      a1 = a2
      u1 =u2
      a2 = 2 * a2
      u2 = fun_dif_mu_zeros(a2,x,somma_zeri)
    }
  }
  a = (a1+a2)/2
  u = fun_dif_mu_zeros(a,x,somma_zeri)
  while (abs(u)>0.001) {
    if(u>0){
      a2 = a
      u2 = u
    }else{
      a1 = a
      u1 = u
    }
    a = (a1 + a2)/2
    u = fun_dif_mu_zeros(a,x,somma_zeri)
  }
  r = data.frame(a,u)
  rownames(r) = x
  return(r)
}


cells_keep <- function(ra,hk, n){
  col_keep = vector()
  col_drop = vector()
  for (i in 1:length(colnames(ra))) {
    if (length(which(ra[hk,i] == 0)) < n ){ 
      col_keep = c(col_keep,i)
    }else{
      col_drop = c(col_drop,i)
    }
  }
  ra = ra[,col_keep]
  return(ra)
}


#### ----- Function to clean data ####

fun_sqrt <- function(cells, raw){    
  # From here the code need to be iterated 
 
  
  print("Start estimation mu")
  # Estimators computation
  raw = raw[rownames(cells),colnames(cells)]
  
  new_raw = sqrt(raw)
 
    print("lambda i")
    lambda_i = lambda(rowMeans(new_raw), new_raw)
    print("nu est")
    nu_est = nu(colMeans(new_raw),new_raw)
    print("mu estimator")
    mu_estimator = outer(lambda_i,nu_est, "*") 
    
   
  gc()
  # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
  to_clust <- as.matrix(raw) %b/% matrix(nu_est, nrow = 1)#raw counts divided for cell efficiency
  colnames(to_clust) = colnames(cells)
  t_to_clust = t(to_clust)
  start_time <- Sys.time()
  
  pca_cells = python_PCA(t_to_clust)
  rownames(pca_cells)=rownames(t_to_clust)  
  end_time <- Sys.time()
  print(paste("pca; time",end_time - start_time, sep = " " ))
  
  #---- Mhalanobis distance
  ppp = pca_cells
  ppp = scale(ppp)
  dist_cells = dist(ppp, method = "euclidean") # mhalanobis
  colnames(pca_cells) = paste("PC",c(1:ncol(pca_cells)), sep = "") 
  pca_cells = as.data.frame(pca_cells)
  output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "nu_est"=nu_est, "lambda_i"=lambda_i)#, "si_si"=si_si, "si_no"=si_no,"no_no"=no_no,"no_si"=no_si,"genes_on_off"=genes_on_off )
  
  return(output)
}

fun_sqrt_iter <- function(cells, raw){    
  print("Start estimation mu sqrt iterative")
  # Estimators computation
  raw = raw[rownames(cells),colnames(cells)]
  
  new_raw = sqrt(raw)
  genes_on_off = as.data.frame(matrix(ncol = length(colnames(new_raw)), nrow = length(rownames(new_raw)),data = TRUE))
  rownames(genes_on_off) = rownames(new_raw)
  colnames(genes_on_off) = colnames(new_raw)
  
  free_deg = 0
  old_free_deg = 1
  while (free_deg != old_free_deg) {
   old_free_deg = free_deg
  print("lambda i")
  lambda_i = lambda(rowMeans(new_raw), new_raw)
  print("nu est")
  nu_est = nu(colMeans(new_raw),new_raw)
  print("mu estimator")
  mu_estimator = outer(lambda_i,nu_est, "*") 
  
   genes_on_off = genes_on_off & !(exp(-mu_estimator) < 10**-3 & raw[rownames(cells),colnames(cells)]  == 0 )
      free_deg = sum(genes_on_off == FALSE, na.rm = T)
     print(free_deg)
  new_raw[!genes_on_off] = sqrt(mu_estimator[!genes_on_off])
     gc()
  
    }
  gc()
  # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
  to_clust <- as.matrix(raw) %b/% matrix(nu_est, nrow = 1)#raw counts divided for cell efficiency
  colnames(to_clust) = colnames(cells)
  t_to_clust = t(to_clust)
  start_time <- Sys.time()
  
  
  pca_cells = python_PCA(t_to_clust)
  rownames(pca_cells)=rownames(t_to_clust)  
  end_time <- Sys.time()
  print(paste("pca; time",end_time - start_time, sep = " " ))
  
  #---- Mhalanobis distance
  ppp = pca_cells
  ppp = scale(ppp)
  dist_cells = dist(ppp, method = "euclidean") # mhalanobis
  colnames(pca_cells) = paste("PC",c(1:ncol(pca_cells)), sep = "") 
  pca_cells = as.data.frame(pca_cells)
  output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "nu_est"=nu_est, "lambda_i"=lambda_i)#, "si_si"=si_si, "si_no"=si_no,"no_no"=no_no,"no_si"=no_si,"genes_on_off"=genes_on_off )
  
  return(output)
  }


fun_linear <- function(cells, raw){    
  print("Start estimation mu with average method")
  raw = as.matrix(raw[rownames(cells),colnames(cells)])
  
  start_time <- Sys.time()
  
    cells_means = colMeans(raw, dims = 1, na.rm = T)
    genes_means = rowMeans(raw, dims = 1, na.rm = T)
    means = mean(raw,na.rm = T )
    
    mu_estimator = (genes_means %*% t(cells_means)) /means
    
    rownames(mu_estimator) = names(genes_means)
    
    nu_est = colMeans(mu_estimator)/means
    lambda_i = genes_means
    
  end_time <- Sys.time()
  print(paste("End estimation; time",end_time - start_time, sep = " " ))
  gc()
  # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
  start_time <- Sys.time()
  
  to_clust <- raw %b/% matrix(nu_est, nrow = 1)
  
  end_time <- Sys.time()
  print(paste("to clust; time",end_time - start_time, sep = " " ))
  
  t_to_clust = t(to_clust)
  
    
    start_time <- Sys.time()
    
  
  pca_cells = python_PCA(t_to_clust)
  rownames(pca_cells)=rownames(t_to_clust)  
  end_time <- Sys.time()
  print(paste("pca; time",end_time - start_time, sep = " " ))
  #dist_cells = dist(t_to_clust, method = "euclidean")
  #---- Mhalanobis distance
  ppp = pca_cells
  ppp = scale(ppp)
  dist_cells = dist(ppp, method = "euclidean") # mhalanobis
  colnames(pca_cells) = paste("PC",c(1:ncol(pca_cells)), sep = "") 
  pca_cells = as.data.frame(pca_cells)
  output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "nu_est"=nu_est, "lambda_i"=lambda_i)#, "si_si"=si_si, "si_no"=si_no,"no_no"=no_no,"no_si"=no_si,"genes_on_off"=genes_on_off )
  return(output)
}


fun_linear_iter <- function(cells, raw){    
  print("Start estimation mu with average method")
  
  raw = as.matrix(raw[rownames(cells),colnames(cells)])
  new_raw = raw
  
  genes_on_off = (matrix(ncol = length(colnames(new_raw)), nrow = length(rownames(new_raw)),data = TRUE))
  rownames(genes_on_off) = rownames(new_raw)
  colnames(genes_on_off) = colnames(new_raw)
  lambda_i = rowMeans(new_raw, dims = 1, na.rm = T)
  free_deg = 0
  old_free_deg = 1
  start_time <- Sys.time()
  it =1
  while (free_deg != old_free_deg) {
    start_time_it <- Sys.time()  
    old_free_deg = free_deg
    
    cells_means = colMeans(new_raw, dims = 1, na.rm = T)
    genes_means = rowMeans(new_raw, dims = 1, na.rm = T)
    means = mean(new_raw,na.rm = T )
    
  end_time <- Sys.time()
  print(paste("after means evaluation",end_time - start_time, sep = " "))
  mu_estimator = (genes_means %*% t(cells_means)) /means
  end_time <- Sys.time()
  print(paste("after putting together",end_time - start_time, sep = " "))
  rownames(mu_estimator) = names(genes_means)
  colnames(mu_estimator) = names(cells_means)
  #print(mu_estimator[1:3,1:3])
  mu_estimator = as.data.frame(mu_estimator)
  #print(mu_estimator[1:3,1:3])
   genes_on_off = genes_on_off & !(exp(-mu_estimator) < 10**-3 & raw[rownames(cells),colnames(cells)] == 0 ) #& raw[rownames(cells),] == 0 )
  
    free_deg = sum(genes_on_off == FALSE, na.rm = T)
   print(paste("free deg",free_deg, sep = " "))
  
    new_raw[!genes_on_off] = mu_estimator[!genes_on_off]
  nu_est = colMeans(mu_estimator)/means
  
  end_time <- Sys.time()
    print(paste("after iteration",it,"time",end_time - start_time_it, sep = " "))
    it = it + 1
    gc()
  }
  end_time <- Sys.time()
  print(paste("End estimation; time",end_time - start_time, sep = " " ))
  gc()
  # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
  start_time <- Sys.time()
  
  to_clust <- raw %b/% matrix(nu_est, nrow = 1)
 
  end_time <- Sys.time()
  print(paste("to clust; time",end_time - start_time, sep = " " ))
  
  t_to_clust = t(to_clust)
  
  start_time <- Sys.time()
  
  
  pca_cells = python_PCA(t_to_clust)
  rownames(pca_cells)=rownames(t_to_clust)  
  end_time <- Sys.time()
  print(paste("pca; time",end_time - start_time, sep = " " ))
  #---- Mhalanobis distance
  ppp = pca_cells
  ppp = scale(ppp)
  dist_cells = dist(ppp, method = "euclidean") # mhalanobis
  colnames(pca_cells) = paste("PC",c(1:ncol(pca_cells)), sep = "") 
  pca_cells = as.data.frame(pca_cells)
  output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "nu_est"=nu_est, "lambda_i"=lambda_i)#, "si_si"=si_si, "si_no"=si_no,"no_no"=no_no,"no_si"=no_si,"genes_on_off"=genes_on_off )
  return(output)
}


fun_linear_iter2 <- function(cells, raw, out_dir){    
  print("Start estimation mu with average method")
  
  raw = as.matrix(raw[rownames(cells),colnames(cells)])
  new_raw = raw
  
  genes_on_off = (matrix(ncol = length(colnames(new_raw)), nrow = length(rownames(new_raw)),data = TRUE))
  rownames(genes_on_off) = rownames(new_raw)
  colnames(genes_on_off) = colnames(new_raw)
  
  free_deg = 0
  old_free_deg = 1
  start_time <- Sys.time()
  it =1
  while (free_deg != old_free_deg) {
    start_time_it <- Sys.time()  
    old_free_deg = free_deg
    
    cells_means = colMeans(new_raw, dims = 1, na.rm = T)
    genes_means = rowMeans(new_raw, dims = 1, na.rm = T)
    means = mean(new_raw,na.rm = T )
    
    end_time <- Sys.time()
    print(paste("after means evaluation",end_time - start_time, sep = " "))
    mu_estimator = (genes_means %*% t(cells_means)) /means
    end_time <- Sys.time()
    print(paste("after putting together",end_time - start_time, sep = " "))
    rownames(mu_estimator) = names(genes_means)
    colnames(mu_estimator) = names(cells_means)
    #print(mu_estimator[1:3,1:3])
    mu_estimator = as.data.frame(mu_estimator)
    #print(mu_estimator[1:3,1:3])
    genes_on_off = genes_on_off & !(exp(-mu_estimator) < 10**-3 & raw[rownames(cells),colnames(cells)] == 0 ) #& raw[rownames(cells),] == 0 )
    
    free_deg = sum(genes_on_off == FALSE, na.rm = T)
    print(paste("free deg",free_deg, sep = " "))
    
    new_raw[!genes_on_off] = mu_estimator[!genes_on_off]
    nu_est = colMeans(mu_estimator)/means
    lambda_i = genes_means
    end_time <- Sys.time()
    print(paste("after iteration",it,"time",end_time - start_time_it, sep = " "))
    it = it + 1
    gc()
  }
  end_time <- Sys.time()
  write.csv(new_raw,paste(out_dir,"new_raw.csv", sep = ""))
  print(paste("End estimation; time",end_time - start_time, sep = " " ))
  gc()
  # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
  start_time <- Sys.time()
  
  to_clust <- raw %b/% matrix(nu_est, nrow = 1)
  
  end_time <- Sys.time()
  print(paste("to clust; time",end_time - start_time, sep = " " ))
 
  t_to_clust = t(to_clust)
  
  start_time <- Sys.time()
  
  
  pca_cells = python_PCA(t_to_clust)
  rownames(pca_cells)=rownames(t_to_clust)  
  end_time <- Sys.time()
  print(paste("pca; time",end_time - start_time, sep = " " ))
  #---- Mhalanobis distance
  ppp = pca_cells
  ppp = scale(ppp)
  dist_cells = dist(ppp, method = "euclidean") # mhalanobis
  colnames(pca_cells) = paste("PC",c(1:ncol(pca_cells)), sep = "") 
  pca_cells = as.data.frame(pca_cells)
  output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "nu_est"=nu_est, "lambda_i"=lambda_i)#, "si_si"=si_si, "si_no"=si_no,"no_no"=no_no,"no_si"=no_si,"genes_on_off"=genes_on_off )
  return(output)
}


fun_iter_no_nu <- function(cells, raw){    
  si_si = as.matrix(cells) %*% t(as.matrix(cells))
  si_si = as.data.frame(si_si)
  
  somma = rowSums(cells)
  somma = as.data.frame(somma)
  si_any = do.call("cbind", replicate(length(rownames(somma)), somma, simplify = FALSE))
  rm(somma)
  colnames(si_any) = rownames(si_any)
  si_no = si_any - si_si
  si_any = t(si_any)
  si_any = as.data.frame(si_any)
  no_si = si_any - si_si
  rm(si_any)
  no_no = n_cells - (si_si + no_si + si_no)
  
  print("Start estimation mu")
  # Estimators computation
  raw = raw[rownames(cells),colnames(cells)]
  
  new_raw = sqrt(raw)
  genes_on_off = as.data.frame(matrix(ncol = length(colnames(new_raw)), nrow = length(rownames(new_raw)),data = TRUE))
  rownames(genes_on_off) = rownames(new_raw)
  colnames(genes_on_off) = colnames(new_raw)
  
  free_deg = 0
  old_free_deg = 1
  while (free_deg != old_free_deg) {
    old_free_deg = free_deg
    
    lambda_i = lambda(rowMeans(new_raw), new_raw)
    nu_est = nu(colMeans(new_raw),new_raw)
    
    mu_estimator = outer(lambda_i,nu_est, "*") 
    
    genes_on_off = genes_on_off & !(exp(-mu_estimator) < 10**-3 & raw[rownames(cells),colnames(cells)]  == 0 )
    free_deg = sum(genes_on_off == FALSE, na.rm = T)
    print(free_deg)
    new_raw[!genes_on_off] = sqrt(mu_estimator[!genes_on_off])
    
  }
  # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
  to_clust = as.matrix(raw[,colnames(cells)])
  colnames(to_clust) = colnames(cells)
  t_to_clust = t(to_clust)
  pca_cells = prcomp(t_to_clust)
  
  #---- Mhalanobis distance
  ppp = pca_cells$x[,1:20] # I try 20 component instead 10
  ppp = scale(ppp)
  dist_cells = dist(ppp, method = "euclidean") # mhalanobis
  
  output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "nu_est"=nu_est, "si_si"=si_si, "si_no"=si_no,"no_no"=no_no,"no_si"=no_si,"genes_on_off"=genes_on_off )
  return(output)
}

#################

