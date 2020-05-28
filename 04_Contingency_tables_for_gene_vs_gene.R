#########################################
# Contingency tables data production for gene vs gene analysis
#########################################
# To run after the data cleaning

########
# Part 1
########

#out_dir = "../results/2019.12.05/"
# ------------------------------------------
{si_si = as.matrix(cells) %*% t(as.matrix(cells))
somma = rowSums(cells)
somma = as.matrix(somma)
si_any = do.call("cbind", replicate(length(rownames(somma)), somma, simplify = FALSE))
rm(somma)
colnames(si_any) = rownames(si_any)
si_no = si_any - si_si
si_any = t(si_any)
no_si = si_any - si_si
rm(si_any)
no_no = n_cells - (si_si + no_si + si_no)

print("Writing contingency tables for observed data.")

#fwrite(as.data.frame(si_si), paste(out_dir,"contingency_table_si_si_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
#fwrite(as.data.frame(si_no), paste(out_dir,"contingency_table_si_no_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
#fwrite(as.data.frame(no_si), paste(out_dir,"contingency_table_no_si_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
#fwrite(as.data.frame(no_no), paste(out_dir,"contingency_table_no_no_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
}

#---------------------------------------------------
####################################
#     ONLY if the observed tables are not already loaded and if they were saved
###################################
cells = as.data.frame(fread(paste(out_dir,"cells_",t,".csv", sep = "")))
rownames(cells)=cells$V1
cells = as.matrix(cells[,2:ncol(cells)])

load(paste(out_dir,"mu_estimator",t, sep = ""))
print("Reading reads tables ")
si_si = as.data.frame(fread(paste(out_dir,"contingency_table_si_si_",t,".csv.gz", sep = "")))
rownames(si_si)=si_si[,1]
si_si = as.matrix(si_si[,2:ncol(si_si)])

si_no =as.data.frame(fread(paste(out_dir,"contingency_table_si_no_",t,".csv.gz", sep = "")))
rownames(si_no)=si_no[,1]
si_no = as.matrix(si_no[,2:ncol(si_no)])

no_si = as.data.frame(fread(paste(out_dir,"contingency_table_no_si_",t,".csv.gz", sep = "")))
rownames(no_si)=no_si[,1]
no_si = as.matrix(no_si[,2:ncol(no_si)])

no_no = as.data.frame(fread(paste(out_dir,"contingency_table_no_no_",t,".csv.gz", sep = "")))
rownames(no_no)=no_no[,1]
no_no = as.matrix(no_no[,2:ncol(no_no)])
#----

########
# Part 2
########

# to run always
{ hk = names(which(rowSums(cells) == length(colnames(cells))))

# exlude the effective ubiqutarius genes and saved in a separate file
mu_estimator = mu_estimator[!rownames(mu_estimator) %in% hk,]
cells = cells[!rownames(cells) %in% hk, ]
write.csv(hk, paste(out_dir,"housekeeping_genes_",t,".csv", sep = ""))

mu_estimator = as.matrix(mu_estimator)
print(paste("start a minimization: time", t, sep = " "))
gc()
p=1
tot = list()
while(p <= length(rownames(mu_estimator))) {
  bb=Sys.time()
  if((p+200) <= length(rownames(mu_estimator))){
    tot1 =  mclapply(rownames(mu_estimator)[p:(p+200)], fun_my_opt,mc.cores = 11)
    
  }else{
    print("Final trance!")
    tot1 = mclapply(rownames(mu_estimator)[p:length(rownames(mu_estimator))], fun_my_opt,mc.cores = 11)
  }
  tot = append(tot, tot1)
  p=p+200+1
  print(Sys.time()-bb)
  print(paste("Next gene:",rownames(mu_estimator)[p],"number",p, sep = " "))
}
gc()
 
tot2 = tot[[1]]
for (tt in 2:length(tot)) {
  tot2= rbind(tot2,tot[[tt]])
  
}

save(tot2 , file = paste(out_dir,"a_minimization_",t, sep = ""))
print(paste("a min:",min(tot2$a) ,"| a max",max(tot2$a) , "| negative a %:",sum(tot2$a <0)/nrow(tot2)*100 ,sep=" "))


print("Matrixes estimator are going to be created!")

M = fun_pzero(tot2$a,mu_estimator[,colnames(cells)])
N = 1-M #fun_pzero(as.numeric(tot2[,2]),mu_estimator[,colnames(cells)])

n_zero_esti = rowSums(M) # estimated number of zeros for each genes
n_zero_obs = rowSums(cells == 0) # observed number of zeros for each genes
dist_zeros = sqrt(sum((n_zero_esti - n_zero_obs)^2)) 

print(paste("The distance between estimated n of zeros and observed number of zero is", dist_zeros,"over", length(rownames(M)), sep = " "))

if(any(is.na(M))){
  print(paste("Errore: some Na in matrix M", which(is.na(M),arr.ind = T),sep = " "))
  break()
}

rm(tot2, mu_estimator)
gc()
estimator_no_no = M %*% t(M)
estimator_no_si = M %*% t(N)
estimator_si_no = t(estimator_no_si)
estimator_si_si = N %*% t(N)

print("Done")
rm(M)
rm(N)

print("Starting dif calculations")

no_si= no_si[!rownames(no_si) %in% hk,!colnames(no_si) %in% hk]
no_no = no_no[!rownames(no_no) %in% hk,!colnames(no_no) %in% hk]
si_si = si_si[!rownames(si_si) %in% hk,!colnames(si_si) %in% hk]
si_no = si_no[!rownames(si_no) %in% hk,!colnames(si_no) %in% hk]

new_estimator_no_no = estimator_no_no
new_estimator_no_no[new_estimator_no_no < 1] <- 1
new_estimator_no_si = estimator_no_si
new_estimator_no_si[new_estimator_no_si < 1] <- 1
new_estimator_si_si = estimator_si_si
new_estimator_si_si[new_estimator_si_si < 1] <- 1
new_estimator_si_no = estimator_si_no
new_estimator_si_no[new_estimator_si_no < 1] <- 1

dif_no_no = (as.matrix(no_no) - estimator_no_no)**2/new_estimator_no_no
dif_no_si = (as.matrix(no_si) - estimator_no_si)**2/new_estimator_no_si
dif_si_no = (as.matrix(si_no) - estimator_si_no)**2/new_estimator_si_no
dif_si_si = (as.matrix(si_si) - estimator_si_si)**2/new_estimator_si_si

S = dif_no_si+dif_si_no+dif_si_si+dif_no_no

if(any(is.na(S))){
  print(paste("Errore: some Na in matrix S", which(is.na(S),arr.ind = T),sep = " "))
  break()
}

rm(dif_no_no)
rm(dif_si_si)
rm(dif_no_si)
rm(dif_si_no)
gc()

fwrite(as.data.frame(S), paste(out_dir,"S_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
print("Calculating p values")

p_value = pchisq(as.matrix(S), df=1, lower.tail=F)
fwrite(as.data.frame(p_value), paste(out_dir,"p_value_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")

rm(p_value)
rm(S)
gc()

new_estimator_si_si = as.matrix(estimator_si_si)
new_estimator_si_si[new_estimator_si_si < 1] <- 1
new_estimator_si_no = as.matrix(estimator_si_no)
new_estimator_si_no[new_estimator_si_no < 1] <- 1
new_estimator_no_no = as.matrix(estimator_no_no)
new_estimator_no_no[new_estimator_no_no < 1] <- 1
new_estimator_no_si = as.matrix(estimator_no_si)
new_estimator_no_si[new_estimator_no_si < 1] <- 1

coex = ((as.matrix(si_si) - as.matrix(estimator_si_si))/new_estimator_si_si) + 
  ((as.matrix(no_no) - as.matrix(estimator_no_no))/new_estimator_no_no) - 
  ((as.matrix(si_no) - as.matrix(estimator_si_no))/new_estimator_si_no) - 
  ((as.matrix(no_si) - as.matrix(estimator_no_si))/new_estimator_no_si)

coex = coex / sqrt(1/new_estimator_si_si + 1/new_estimator_no_no + 1/new_estimator_si_no + 1/new_estimator_no_si)

fwrite(as.data.frame(coex), paste(out_dir,"coex_",t,".csv.gz", sep = ""),row.names = T)
rm(coex)
gc()


#fwrite(as.data.frame(estimator_si_si), paste(out_dir,"estimator_si_si_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
#fwrite(as.data.frame(estimator_si_no), paste(out_dir,"estimator_si_no_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
#fwrite(as.data.frame(estimator_no_si), paste(out_dir,"estimator_no_si_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
#fwrite(as.data.frame(estimator_no_no), paste(out_dir,"estimator_no_no_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")

rm(estimator_si_si)
rm(estimator_si_no)  
rm(estimator_no_si)  
rm(estimator_no_no)  

rm(new_estimator_si_si)
rm(new_estimator_si_no)  
rm(new_estimator_no_si)  
rm(new_estimator_no_no)


rm(si_si)
rm(si_no)
rm(no_si)
rm(no_no)


gc()
}
