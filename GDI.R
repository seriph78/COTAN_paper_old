source("../src/functions.R")
library(ggrepel)

######

### Global Differentiation Index 
# Set some genes
hk2 = c("Calm1","Cox6b1","Ppia","Rpl18","Cox7c","Erh","H3f3a","Tbp","Taf1","Taf2","Gapdh","Actb")
tf1 = c("Nes","Vim","Sox2","Sox1","Notch1", "Hes1","Hes5","Hes3","Pax6")
tf2 = c("Map2","Tubb3","Neurod1","Nefm","Nefl","Dcx","Tbr1")
layer = c("Reln","Lhx5","Cux1","Satb2","Tle1","Mef2c","Rorb","Sox5","Bcl11b","Fezf2","Foxp2")

genes = c("Gapdh","Actb","Tbp","Satb2","Vim","Hes1","Bcl11b")
######

# use the function that calculates the GDI: it needs the t prefix for the dataset, the matrixes folder
# and the code for the group of dataset 

quant.p = GDI_dataframe(dir = "../results/2020.01.26/",t = 17,s = "M1")

quant.p.sorted = quant.p[order(quant.p$GDI, decreasing = T),]


quant.p$highlight = with(quant.p, ifelse(rownames(quant.p) %in% tf1, "tf1",
                                                 ifelse(rownames(quant.p) %in% hk2,"hk" ,
                                                        ifelse(rownames(quant.p) %in% layer,"layer" ,
                                                               ifelse(rownames(quant.p) %in% tf2,"tf2" , "normal"))))) 
textdf <- quant.p[rownames(quant.p) %in% c(tf1,hk2,layer,tf2), ]
mycolours <- c("hk" = "#3C5488B2","tf1"="#F39B7FE5","layer"="#7E6148B2","tf2"="#E64B35B2","normal" = "#91D1C2B2")


#pdf(file = paste(out_dir,"GDI_iteration_new",t, ".pdf", sep = ""), paper = "a4r")
ggplot(quant.p, aes(x = sum.raw.norm, y = GDI)) +
  geom_point(size = 2, aes(colour = highlight), alpha=0.5) +
  scale_color_manual("Status", values = mycolours) +
#  ylim(0,20.5)+
  geom_label_repel(data = textdf, aes(x = sum.raw.norm, y = GDI, label = rownames(textdf)),hjust=0.5, vjust=-0.5) +
  geom_hline(yintercept=1.5)+#, linetype="dashed", color = "red", size 2) +
  geom_hline(yintercept=quantile(quant.p$GDI)[4], linetype="dashed", color = "darkblue") +
  geom_hline(yintercept=quantile(quant.p$GDI)[3], linetype="dashed", color = "darkblue") +
  #geom_vline(xintercept=quantile(sum.coex$exp.introns)[4], linetype="dashed", color = "darkred") +
  #geom_vline(xintercept=quantile(sum.coex$exp.introns)[3], linetype="dashed", color = "darkred") +
  xlab("log normalized reads sum")+
  ylab("global p val index (GDI)")+
  theme(legend.position = "none") +
#  ggtitle(paste("GDI ",t, sep = ""))+
  theme()
#dev.off()

#fwrite(quant.p.val, paste(out_dir,"quant_p_val_",t,".csv.gz", sep = ""),row.names = T,compress = "gzip")
