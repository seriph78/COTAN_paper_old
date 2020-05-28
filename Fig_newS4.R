# Figure new S4
source("../src/functions.R")
source("../src/cotan_output_functions.R")
library("ggsci")
library("ggplot2")
data.dir = "../results/2020.01.26/"
primary.markers = c("Satb2","Bcl11b","Sox5","Reln","Rorb","Tbr1","Lhx5","Cux1","Fezf2","Foxp2","Vim","Hes1")

load_data3.0(dir = data.dir, cond = 17,genes = primary.markers ,prefix = "p_value_" )

n.genes.for.marker = 25
all.genes.to.an = vector()
for (m in primary.markers) {
  all.genes.to.an = c(all.genes.to.an,rownames(p_value_17[order(p_value_17[,m]),])[1:n.genes.for.marker])
  all.genes.to.an =unique(all.genes.to.an)
}

load_data3.0(dir = data.dir, cond = 17 ,genes = all.genes.to.an ,prefix = "p_value_" )

S =as.data.frame(fread(paste(data.dir,"S_17.csv.gz", sep="")) )
#S = cbind(S,as.data.frame(fread(paste(data.dir,"S_17.csv.gz", sep=""),select = all.genes.to.an)))
rownames(S)= S$V1
S=S[,2:ncol(S)]

#quant.p.val2 = rowQuantiles((as.matrix(S)),probs =(1-15/ncol(S)) , na.rm = T) #0.975
CD.sorted <- t(apply(t(S[,colnames(S) %in% all.genes.to.an]),2,sort,decreasing=T))
CD.sorted = CD.sorted[,1:round(ncol(CD.sorted)/20, digits = 0)]
CD.sorted = pchisq(as.matrix(CD.sorted), df=1, lower.tail=F)

GDIpar = rowMeans(CD.sorted)
GDIpar =as.data.frame(GDIpar)
colnames(GDIpar) = "S"


#GDIpar =as.data.frame(GDIpar)
#colnames(GDIpar) = "S"
GDIpar$names = rownames(GDIpar)
GDIpar$log.mean.p = -log(GDIpar$S)
GDIpar$GDIpar = log(GDIpar$log.mean.p)

#########################
# GDI GW
CD.sorted <- apply(S,2,sort,decreasing=T)
CD.sorted = CD.sorted[1:round(nrow(CD.sorted)/20, digits = 0),]
CD.sorted = pchisq(as.matrix(CD.sorted), df=1, lower.tail=F)

GDIgw = colMeans(CD.sorted)
GDIgw =as.data.frame(GDIgw)
colnames(GDIgw) = "S"

GDIgw$log.mean.p = -log(GDIgw$S)
GDIgw$GDIgw = log(GDIgw$log.mean.p)

GDI = merge(GDIgw, GDIpar, by="row.names", all=TRUE) 


# just to color
layers = list("L1"=c("Reln","Lhx5"), "L2"=c("Satb2","Cux1"), "L4"=c("Rorb","Sox5") , "L5"=c("Bcl11b","Fezf2") , "P"=c("Vim", "Hes1"))

pos.link = list()
for (m in c(1:length(layers))) {
  print(m)
  tmp2 = vector()
  for (g in layers[[m]]) {
    print(g)
    tmp =rownames(p_value_17[order(p_value_17[,g]),])[1:n.genes.for.marker+20]
    tmp2= c(tmp2,tmp[tmp %in% rownames(coex_17[coex_17[,g] >= 0,])]  )
  }
  pos.link[[m]] = unique(c(tmp2,layers[[m]]))
}
names(pos.link)=names(layers)
dup = unlist(pos.link)[duplicated(unlist(pos.link))]
# remove the duplicated genes
for (m in c(1:length(pos.link))) {
  pos.link[[m]]= pos.link[[m]][!pos.link[[m]] %in% dup]
  if (!all(layers[[m]] %in% pos.link[[m]])) {
    layers[[m]][!layers[[m]] %in% pos.link[[m]]]
  }
  
}



#genes.raw = S[S >= quant.p.val2$S,]$names #0.9
#genes.raw[grep('^(X[0-9])', genes.raw)] = substring(genes.raw[grep('^(X[0-9])', genes.raw)], 2)
#genes[grep("-",genes)] =  gsub("-",".",genes[grep("-",genes)])
rownames(GDI) = GDI$Row.names

GDI$highlight = with(GDI, 
                       ifelse(rownames(GDI) %in% pos.link$L5, "genes related to L5/6",
                          ifelse(rownames(GDI) %in% pos.link$L2 , "genes related to L2/3",
                          ifelse(rownames(GDI) %in% pos.link$P , "genes related to Prog" ,
                          ifelse(rownames(GDI) %in% pos.link$L1 , "genes related to L1" ,
                          ifelse(rownames(GDI) %in% pos.link$L4 ,"genes related to L4" ,
                                                          "not marked"))))))

controls =list("genes related to L5/6"=c("Tle4","Foxp2","Tbr1"), "genes related to L2/3"=c("Mef2c"), "genes related to Prog"=c("Nes","Sox2") , "genes related to L1"=c() , "genes related to L4"=c("Kcnip2")) 


textdf <- GDI[rownames(GDI) %in% c(unlist(layers),unlist(controls)) , ]

for (m in c(1:length(controls))) {
  for (g in controls[[m]]) {
    if(g %in% rownames(textdf)){
      textdf[g,"highlight"] = names(controls[m])
    } 
  }
}



mycolours <- c("genes related to L5/6" = "#3C5488FF","genes related to L2/3"="#F39B7FFF","genes related to Prog"="#4DBBD5FF","genes related to L1"="#E64B35FF","genes related to L4" = "#91D1C2FF", "not marked"="#B09C85FF")

tmp1 = ggplot(subset(GDI,highlight == "not marked" ), aes(y=GDIgw, x=GDIpar)) +  geom_point(alpha = 0.3, color = "#B09C85FF",size=2.5)

GDIplot = tmp1 + geom_point(data = subset(GDI, highlight != "not marked" ), aes(y=GDIgw, x=GDIpar, colour=highlight),size=2.5,alpha = 0.9) +
  scale_color_manual("Status", values = mycolours)  +
 # xlim(0,5)+
#  ylim(0,5)+
  scale_fill_manual("Status", values = mycolours)  +
  ylab("genome-wide GDI") + xlab("localized GDI")+
  geom_abline(slope = 1,linetype="dotted", colour = "darkred")+ 
  geom_label_repel(data =textdf , aes(y=GDIgw, x=GDIpar, label = rownames(textdf),fill=highlight),
                   label.size = NA, 
                   alpha = 0.5, 
                   direction = "both",
                   na.rm=TRUE,
                   seed = 1234) +
  geom_label_repel(data =textdf , aes(y=GDIgw, x=GDIpar, label = rownames(textdf)),
                   label.size = NA, 
                   segment.color = 'black',
                   segment.size = 0.5,
                   direction = "both",
                   alpha = 1, 
                   na.rm=TRUE,
                   fill = NA,
                   seed = 1234) +
  #ggtitle("PCA") +
  theme_light(base_size=10) +
  theme(axis.text=element_text(size=12, 
                                 face="italic", 
                                 color="#3C5488FF"),
        plot.title = element_text(size=14, 
                                  face="italic", 
                                  color="#3C5488FF",
                                  hjust=0.01,
                                  lineheight=1.2,margin = margin(t = 5, b = -15)),
        axis.title =element_text(size=12, 
                                 face="italic", 
                                 color="#3C5488FF"),
        legend.position = "none")  # titl)


si = 12

themex= theme(axis.text.x = element_text( size = si, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
              axis.text.y = element_blank(),  
              axis.title.x = element_blank(),
              axis.title.y = element_text( size = si, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
              legend.position = "none") 

themey = theme(axis.text.y = element_text( size = si, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
               #axis.text.y = element_blank(),  
               axis.title.x = element_text( size = si, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
               axis.text.x.bottom = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "none") 



xdensityGDI <- ggplot(GDI, aes(GDIpar)) + 
  geom_density(alpha=.5, fill = "#8491B4B2", colour ="#8491B4B2" ) +
  themex  +geom_vline(aes(xintercept=median(GDIpar)),
                      color="#8491B4B2", linetype="dashed", size=1)

ydensityGDI <- ggplot(GDI, aes(GDIgw)) + 
  geom_density(alpha=.5, fill="#00A087FF", colour= "#00A087FF") + 
  themey +  coord_flip()+ geom_vline(aes(xintercept=median(GDIgw)),
                                    color="#00A087FF", linetype="dashed", size=1)


GDIplot = xdensityGDI + plot_spacer() + GDIplot+ ydensityGDI + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))



GDIplot



quant.p = GDI_dataframe(t="17",data.dir,"M1")


quant.p$highlight = with(quant.p, ifelse(rownames(quant.p) %in% pos.link$L5, "genes related to L5/6",
                                         ifelse(rownames(quant.p) %in% pos.link$L2 , "genes related to L2/3",
                                                ifelse(rownames(quant.p) %in% pos.link$P , "genes related to Prog" ,
                                                       ifelse(rownames(quant.p) %in% pos.link$L1 , "genes related to L1" ,
                                                              ifelse(rownames(quant.p) %in% pos.link$L4 ,"genes related to L4" ,
                                                                     "not marked"))))))

textdf <- quant.p[rownames(quant.p) %in% c(unlist(layers),unlist(controls)), ]
#mycolours <- c("Constitutive" = "#00A087FF","NPGs"="#E64B35FF","PNGs"="#F39B7FFF","normal" = "#8491B4B2")
#"LIGs"="#7E6148B2",

f1 = ggplot(subset(quant.p,highlight == "normal" ), aes(x=sum.raw.norm, y=GDI)) +  geom_point(alpha = 0.1, color = "#8491B4B2", size=2.5)

si = 12

GDI_plot2 = f1 +  geom_point(data = subset(quant.p,highlight != "normal"  ), aes(x=sum.raw.norm, y=GDI, colour=highlight),size=2.5,alpha = 0.8) +
  geom_hline(yintercept=quantile(quant.p$GDI)[4], linetype="dashed", color = "darkblue") +
  geom_hline(yintercept=quantile(quant.p$GDI)[3], linetype="dashed", color = "darkblue") +
  geom_hline(yintercept=1.5, linetype="dotted", color = "red", size= 0.5) +
  scale_color_manual("Status", values = mycolours)  +
  scale_fill_manual("Status", values = mycolours)  +
  xlab("log normalized counts")+ylab("GDI")+
  geom_label_repel(data =textdf , aes(x=sum.raw.norm, y=GDI, label = rownames(textdf),fill=highlight),
                   label.size = NA, 
                   alpha = 0.5, 
                   direction ="both",
                   na.rm=TRUE,
                   seed = 1234) +
  geom_label_repel(data =textdf , aes(x=sum.raw.norm, y=GDI, label = rownames(textdf)),
                   label.size = NA, 
                   segment.color = 'black',
                   segment.size = 0.5,
                   direction = "both",
                   alpha = 0.8, 
                   na.rm=TRUE,
                   fill = NA,
                   seed = 1234) +
  theme(axis.text.x = element_text(size = si, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
        axis.text.y = element_text( size = si, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),  
        axis.title.x = element_text( size = si, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
        axis.title.y = element_text( size = si, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
        legend.title = element_blank(),
        legend.text = element_text(color = "#3C5488FF",face ="italic" ),
        legend.position = "bottom")  # titl)




