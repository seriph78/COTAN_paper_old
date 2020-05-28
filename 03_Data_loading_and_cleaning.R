###############################################
# 01 Raw data loading and cleaning
###############################################
source("../rstudio/src/functions.R")
source_python("../rstudio/src/python_PCA.py")
library(ggrepel)
library(latex2exp)
if(!file.exists("../rstudio/storage/cotan/results/mouse_dentate_gyrus/")){
  dir.create(file.path("../rstudio/storage/cotan/results/mouse_dentate_gyrus/"))
}

out_dir = "../rstudio/storage/cotan/results/mouse_dentate_gyrus/linear/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

if(!file.exists("../rstudio/storage/cotan/results/mouse1/")){
  dir.create(file.path("../rstudio/storage/cotan/results/mouse1/"))
}

out_dir = "../rstudio/storage/cotan/results/mouse1/linear/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}


# t is the analide condition or time point
t = "E17.5"

# row reads count loading
colnames(raw)=paste(t,"_",colnames(raw), sep = "")

#genes_to_rem = rownames(raw[grep('^MT-', rownames(raw)),]) #genes to remove : mithocondrial
genes_to_rem = rownames(raw[grep('^mt', rownames(raw)),]) #genes to remove : mithocondrial

raw = raw[!rownames(raw) %in% genes_to_rem,]
cells_to_rem = colnames(raw[which(colSums(raw) == 0)])
raw = raw[,!colnames(raw) %in% cells_to_rem]

print(paste("Condition ",t,sep = ""))
#--------------------------------------
n_cells = length(colnames(raw))
print(paste("n cells", n_cells, sep = " "))
cells =raw     
#---------------------------------------------------
# Cells matrix : formed by row data matrix changed to 0-1 matrix
cells[cells > 0] <- 1
cells[cells <= 0] <- 0
# We want to discard genes having less than 3 not 0 counts over 1000 cells
cells = cells[rowSums(cells)> round((length(colnames(raw))/1000*3), digits = 0),] 

#counter for cleaning iteration
n_it = 1

raw = raw[rownames(cells), colnames(cells)]
#------------------------------------------------------------------
# subdirectory formation
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

if(!file.exists(paste(out_dir,"cleaning", sep = ""))){   
  dir.create(file.path(out_dir, "cleaning"))
}

#################################################################
mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
      axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),  
      axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
      axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))
        # titl)

# iterate this part until all outlier cells are removed
# cells  should be a data.frame not a matrix otherwise the function that minimizes a does not work correctly
{ cells = cells[rowSums(cells)> round((length(colnames(raw))/1000*3), digits = 0),] # 3 per mille delle cellule
  raw = raw[rownames(cells), colnames(cells)]
  #fwrite(cells,paste(out_dir,"cells_",t,".csv", sep = ""))
  #save(n_it,file = paste(out_dir,"n_it_",t, sep = ""))
  list1 = fun_linear(cells = cells, raw = raw)
  # list1 = fun_iter_no_nu(cells = cells, raw = raw)
  gc()    
  dist_cells = list1$dist_cells
  pca_cells = list1$pca_cells
  t_to_clust = as.matrix(list1$t_to_clust)
  mu_estimator = list1$mu_estimator
  nu_est = list1$nu_est
  to_clust = list1$to_clust
  
  hc_cells = hclust(dist_cells, method = "complete")
  
  groups <- cutree(hc_cells, k=2)
  
  groups[groups == 1] <- "A"
  groups[groups == 2] <- "B"
  
  cl2 = names(which(cutree(hc_cells, k = 2) == 2))
  cl1 = names(which(cutree(hc_cells, k = 2) == 1))
  
  t_to_clust = cbind(as.data.frame(t_to_clust),groups)
  
  # ---- next: to check which genes are specific for the B group of cells
  B = as.data.frame(to_clust[,colnames(to_clust) %in% cl2])
  colnames(B)=cl2
  B = rownames_to_column(B)
  if (dim(B)[2]>2) {
    B = B[order(rowMeans(B[,2:length(colnames(B))]),decreasing = T), ]
  }else{
    B = B[order(B[,2],decreasing = T), ]  
  }
  
  #B = B[order(B[,2:length(colnames(B))],decreasing = T), ] #if just one column
  
  print(head(B, 15))
  
  C = arrange(B,rowMeans(B[2:length(colnames(B))]))
  rownames(C) = C$rowname
  D = data.frame("means"=rowMeans(C[2:length(colnames(C))]),"n"=NA ) 
  D = D[D$means>0,]
  D$n = c(1:length(D$means))
  
# check if the pca plot is clean enought and from the printed genes, if the smalest group of cells are caratterised by particular genes 

pca_cells = cbind(pca_cells,"groups"=t_to_clust$groups)

#autoplot(pca_cells, data = t_to_clust, colour = 'groups')

pca.cell.1 = ggplot(subset(pca_cells,groups == "A" ), aes(x=PC1, y=PC2,colour =groups)) +geom_point(alpha = 0.5, size=3) 
#  scale_color_manual("groups", values = mycolours)  


pca.cell.2=  pca.cell.1 + geom_point(data = subset(pca_cells, groups != "A" ), aes(x=PC1, y=PC2,colour =groups),alpha = 0.8, size=3)+
  scale_color_manual("groups", values = mycolours)  +
   my_theme + theme(legend.title = element_blank(),
                   legend.text = element_text( size = 12,color = "#3C5488FF",face ="italic" ),
                   legend.position="bottom")
}
pca.cell.2
#---------- run this when B cells are to be removed

pdf(paste(out_dir,"cleaning/",t,"_",n_it,"_plots_before_cells_exlusion.pdf", sep = ""))
pca.cell.2
ggplot(D, aes(x=n,y=means)) + geom_point() + 
  geom_text_repel(data=subset(D, n > (max(D$n)- 15) ), aes(n,means,label=rownames(D[D$n > (max(D$n)- 15),])),
                  nudge_y      = 0.05,
                  nudge_x      = 0.05,
                  direction    = "x",
                  angle        = 90,
                  vjust        = 0,
                  segment.size = 0.2)+ 
  ggtitle("B cell group genes mean expression")+my_theme +
  theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 5,hjust = 0.02 ))
dev.off()

if (length(cl1) < length(cl2)) {
  to_rem = cl1
}else{
  to_rem = cl2
}
n_it = n_it+1
cells = cells[,!colnames(cells) %in% to_rem]
raw = raw[rownames(cells),colnames(cells)]

gc()

#---- run this only in the last iteration, instead the previus code, when B cells group has not to be removed
pdf(paste(out_dir,"cleaning/",t,"_",n_it,"_plots_before_cells_exlusion.pdf", sep = ""))
pca.cell.2
ggplot(D, aes(x=n,y=means)) + geom_point() + 
  geom_text_repel(data=subset(D, n > (max(D$n)- 15) ), aes(n,means,label=rownames(D[D$n > (max(D$n)- 15),])),
                  nudge_y      = 0.05,
                  nudge_x      = 0.05,
                  direction    = "x",
                  angle        = 90,
                  vjust        = 0,
                  segment.size = 0.2)+ 
  ggtitle(label = "B cell group genes mean expression", subtitle = " - B group NOT removed -")+my_theme +
  theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 10,hjust = 0.02 ),
        plot.subtitle = element_text(color = "darkred",vjust = - 15,hjust = 0.01 ))
  
dev.off()
#-----------------------------------------------------------------------

#---------To color the pca based on nu_j (so the cells' efficiency)
nu_est = round(nu_est, digits = 7)

p<-ggplot(pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))

pdf(paste(out_dir,"cleaning/",t,"_plots_PCA_efficiency_colored.pdf", sep = ""))
# or tiff("plot.tiff")
p+geom_point(size = 2.5,alpha= 0.8)+ 
  scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" , 
                        midpoint = mean((nu_est)),name = TeX(" $ln (\\nu) $ "))+
       ggtitle("Cells PCA coloured by cells efficiency") +
  my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                    legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                    legend.text = element_text(color = "#3C5488FF", size = 11),
                    legend.key.width = unit(2, "mm"),
                    legend.position="right")
  

dev.off()

#--------------------
# part to remove the cells with efficiency too low (by hand)
nu_df = data.frame("nu"= sort(nu_est), "n"=c(1:length(nu_est)))

ggplot(nu_df, aes(x = n, y=nu)) + geom_point(colour = "#8491B4B2", size=1)+my_theme + ylim(0,0.7) + xlim(0,250)

pdf(paste(out_dir,"cleaning/",t,"_plots_efficiency.pdf", sep = ""))
# look to the plot and set the treshold where the line bend 
ggplot(nu_df, aes(x = n, y=nu)) + geom_point(colour = "#8491B4B2", size=1)+my_theme + ylim(0,0.7) + xlim(0,250)+
       geom_hline(yintercept=0.4, linetype="dashed", color = "darkred") +
       annotate(geom="text", x=70, y=0.58, label="to remove cells with nu < 0.4 ", color="darkred", size=4.5) 
dev.off()

# or if nothing to remove
pdf(paste(out_dir,"cleaning/",t,"_plots_efficiency.pdf", sep = ""))
ggplot(nu_df, aes(x = n, y=nu)) + geom_point() + ylim(0,1) + xlim(0,500)+ 
       annotate(geom="text", x=50, y=0.25, label="nothing to remove ", color="red") 

dev.off()

# in this case the bend is at 0.25 so the cells having an  estimated efficiency lower than 0.25 are removed
to_rem = rownames(nu_df[which(nu_df$nu < 0.4),]) 

cells = cells[,!colnames(cells) %in% to_rem]
raw = raw[, colnames(cells)]
cells = cells[rowSums(cells)> round((length(colnames(raw))/1000*3), digits = 0),] 
raw = raw[rownames(cells), colnames(cells)]
n_cells = length(colnames(cells))

# repeat the estimation after the cells are removed
list1 = fun_linear(cells = cells, raw = raw)
dist_cells = list1$dist_cells
pca_cells = list1$pca_cells
t_to_clust = list1$t_to_clust   
mu_estimator = list1$mu_estimator
nu_est = list1$nu_est

# to plot the final efficiency colored PCA
p<-ggplot(pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))
# and to save it
pdf(paste(out_dir,"cleaning/",t,"_plots_PCA_efficiency_colored_FINAL.pdf", sep = ""))
p+geom_point(size = 2.5,alpha= 0.8)+ 
  scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" , 
                        midpoint = mean((nu_est)),name = TeX(" $ln (\\nu) $ "))+
  ggtitle("Cells PCA coloured by cells efficiency: last") +
  my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                    legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                    legend.text = element_text(color = "#3C5488FF", size = 11),
                    legend.key.width = unit(2, "mm"),
                    legend.position="right")
dev.off()

# just to check again to keep only genes that are expressed
cells = cells[rowSums(cells)> round((length(colnames(raw))/1000*3), digits = 0),] 
dim(cells)
dim(raw)

save(mu_estimator, file=paste(out_dir,"mu_estimator",t, sep = ""))
fwrite(as.data.frame(cells), paste(out_dir,"cells_",t,".csv", sep = ""),row.names = T)
save(n_cells , file = paste(out_dir,"n_cells",t, sep = ""))
write.csv(nu_est, paste(out_dir,"nu_est_",t,".csv", sep = ""))
write.csv(list1$lambda_i, paste(out_dir,"lambda_i_",t,".csv", sep = ""))

