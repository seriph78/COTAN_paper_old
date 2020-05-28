# function to read the output from COTAN and produce the heatmap
library(data.table)
library(Rfast)
library(dplyr)
library(ggplot2)
library(textshape)

#-------------------------------------
# p_v the desired threshold for p-value
# df_genes is a list of gene sets named. By default the first set will be displayed in the output columns
# sets an array of number defining which of the previous sets ar to be considered
# conditions which are the conditions or time points to plot (small columns)
# ldir directory where are stored the elaborated data (p-values and coex matrix)
# lim-coex: because coex index depends also from the expression and, of course, from the co-expression or
# mutual exclusion, can be useful to set the limits of colour plotting for this index, especially to compare 
# two different plots on the same dataset but with different set of genes (ie. ubiquitous genes and 
# differentiated genes).

heatmap2.0 <- function(p_v, df_genes, sets, conditions, ldir, lim_coex){
  
  
  gr = df_genes[[1]]
  
  ge = array(sort(unlist(df_genes)))
  
  if (0 %in% sets) {
   ge = c()
    for (n in conditions){
    subset_pv = get(paste("p_value_",n, sep = ""))
    if (length(gr) > 1) {
     subset_pv = subset_pv[,colnames(subset_pv) %in% gr]
      subset_pv2 = subset_pv[which(apply(subset_pv, 1, FUN=min) <= p_v),]
      
  }else{
      subset_pv2 = subset_pv[which(subset_pv[,gr] <= p_v),]  
    }
   ge =  unique(c(ge,rownames(subset_pv2)))
    }
  }
  
  to_print_tot = data.frame()
  
  #n = 13
  dl <- data.frame(type = rep(names(df_genes), sapply(df_genes, length)), g2 = unlist(df_genes))
  gg_in_sets <- unique(c(unlist(df_genes[sets]),gr))
  
  #n = "E60"
  #ldir = "results/2019.08.10/"
  
  for (n in conditions){
    #load(paste(ldir,"n_cells",n, sep = ""))
    
    true_hk = read.csv(paste(ldir,"housekeeping_genes_",n,".csv", sep = ""))
    
    true_hk = as.vector(true_hk[true_hk[,2] %in% c(ge,gr),2])
    
    subset_pv = as.matrix(get(paste("p_value_",n, sep = "")))
    
    subset_pv = subset_pv[rownames(subset_pv) %in% ge | rowMins(subset_pv) <= p_v,]
    
    coex = as.matrix(get(paste("coex_",n, sep = "")))
    
    coex = coex[rownames(coex) %in% ge,]
    
    load(paste(ldir,"n_cells",n, sep = ""))
    
    coex = coex/sqrt(n_cells)
    
  #  coex = read.csv(file=paste(ldir,"coex_",n,".csv", sep = ""))
    
    for (b in true_hk){
      v = rep.int(0,times = length(rownames(coex)))
      v = as.data.frame(v)
      colnames(v) =b
      coex = cbind(coex, v)
      
      z = rep.int(0,times = length(colnames(coex)))
      z=as.data.frame(z)
      colnames(z) =b
      rownames(z)=colnames(coex)
      coex = rbind(coex, t(z))
    }
    
    #g1 = ex columns ; g2 = ex rows
    to_print = setNames(reshape2::melt(as.matrix(coex)), c('g2', 'g1', 'coex'))
    to_print$g2 = as.vector(to_print$g2) 
    to_print$g1 = as.vector(to_print$g1) 
    to_print_pv = setNames(reshape2::melt(subset_pv), c('g2', 'g1', 'p_val'))
    to_print_pv$g2 = as.vector(to_print_pv$g2) 
    to_print_pv$g1 = as.vector(to_print_pv$g1) 
    
    #to_print = merge(as.data.frame(to_print), as.data.frame(to_print_pv),  by = c("g1","g2"), all = TRUE)
    to_print <- to_print %>% left_join(to_print_pv, by=c("g1","g2"))
    
    to_print = within(to_print, coex[p_val > p_v] <- 0)
    
    to_print = within(as.data.frame(to_print), coex[as.character(g1) == as.character(g2)] <- 0)
    
    to_print$time = n
    
    if(!0 %in% sets){ # if in sets there isn't the 0, subset only interesting genes to be printed
      to_print = to_print[to_print$g2 %in% ge,]
      to_print = to_print[to_print$g1 %in% gg_in_sets,]
      to_print = merge(to_print, dl, all = T)
      to_print = to_print[to_print$type %in% names(df_genes[sets]),]
    }else{ # 0 means that with consider all genes under p-vales treshold
      to_print = merge(to_print, dl, all = T)
    }
    
    
    to_print = to_print[to_print$g1 %in% df_genes[[1]],]
    
    to_print$t_hk = ifelse((to_print$g2  %in% true_hk) | (to_print$g1  %in% true_hk), "hk", "n")
    
    to_print_tot = rbind(to_print_tot, to_print)
  }
  
  print(paste("min coex:",min(to_print_tot$coex), "max coex", max(to_print_tot$coex),sep = " "))
  
  ggplot(to_print_tot, aes(time, factor(g2, 
                                        levels = rev(levels(factor(g2)))))) + 
    geom_tile(aes(fill = coex),colour = "black", show.legend = TRUE) +
    facet_grid( type ~ g1  ,scales = "free", space = "free") + 
    scale_fill_gradient2(low = "darkred", mid = "white",  high = "darkblue", midpoint = 0,na.value = "grey80", space = "Lab", guide = "colourbar", aesthetics = "fill", limits = lim_coex, oob=scales::squish)+ theme(legend.position="bottom")+
    theme(axis.title.x = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.background = element_rect(fill="#8491B44C"),
          strip.text.y  = element_text(size = 9, colour = "#3C5488FF"),
          strip.text.x  = element_text(size = 9,angle= 90 ,colour = "#3C5488FF"),
          axis.title.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.text.y = element_text( size = 9, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),  
          #legend.title = element_blank(),
          legend.text = element_text(color = "#3C5488FF",face ="italic" ),
          legend.position = "bottom",
          legend.title=element_blank(),
          legend.key.height = unit(2, "mm"))#+geom_text(aes(label=ifelse(t_hk == "hk", "H","")), color="grey", size=3)
  
  
}
#------------------------
# Similar to the previous but the output will be genome wide. 
# The last term (all_conditions) can be True or False
heatmap2.0_gw <- function(p_v,df_genes, genes, conditions, ldir, lim_coex, all_conditions){
  
  gr = genes #vector with the interesting genes
  
  #ge = array(sort(unlist(df_genes)))
  to_print_tot = data.frame()
  
  #n = 13
  dl <- data.frame(type = rep(names(df_genes), sapply(df_genes, length)), g2 = unlist(df_genes))
  gg_in_sets <- gr #unique(c(unlist(df_genes[sets]),gr))
  
    ge = c(unlist(df_genes))
    for (n in conditions){
      subset_pv = get(paste("p_value_",n, sep = ""))
      subset_coex = as.matrix(get(paste("coex_",n, sep = "")))
      
      load(paste(ldir,"n_cells",n, sep = ""))
      
      subset_coex= subset_coex/sqrt(n_cells)
      
      if (length(gr) > 1) {
        subset_pv = subset_pv[,colnames(subset_pv) %in% gr]
        subset_pv = subset_pv[ which(apply(subset_pv, 1, FUN=median) <= p_v | rownames(subset_pv) %in% ge )  ,]
      }else{
        subset_pv = as.data.frame(subset_pv[which(subset_pv[,gr] <= p_v | rownames(subset_pv) %in% ge ),gr, drop=F])
        colnames(subset_pv) = gr
      }
      subset_coex = as.data.frame(subset_coex[rownames(subset_coex) %in% rownames(subset_pv),colnames(subset_coex) %in%  colnames(subset_pv)])
      colnames(subset_coex) =colnames(subset_pv)
      ge =  unique(c(ge,rownames(subset_pv)))
      
      true_hk = read.csv(paste(ldir,"housekeeping_genes_",n,".csv", sep = ""))
      true_hk = as.vector(true_hk[true_hk[,2] %in% c(ge,gr),2])
      
      for (b in true_hk){
        v = rep.int(0,times = length(rownames(coex)))
        v = as.data.frame(v)
        colnames(v) =b
        subset_coex = cbind(subset_coex, v)
        
        z = rep.int(0,times = length(colnames(subset_coex)))
        z=as.data.frame(z)
        colnames(z) =b
        rownames(z)=colnames(subset_coex)
        subset_coex = rbind(subset_coex, t(z))
      }
      to_print = setNames(reshape2::melt(as.matrix(subset_coex)), c('g2', 'g1', 'coex'))
      to_print$g2 = as.vector(to_print$g2) 
      to_print$g1 = as.vector(to_print$g1) 
      to_print_pv = setNames(reshape2::melt(as.matrix(subset_pv)), c('g2', 'g1', 'p_val'))
      to_print_pv$g2 = as.vector(to_print_pv$g2) 
      to_print_pv$g1 = as.vector(to_print_pv$g1) 
      
      #to_print = merge(as.data.frame(to_print), as.data.frame(to_print_pv),  by = c("g1","g2"), all = TRUE)
      to_print <- to_print %>% left_join(to_print_pv, by=c("g1","g2"))
      
      to_print = within(to_print, coex[p_val > p_v] <- 0)
      
      to_print = within(as.data.frame(to_print), coex[as.character(g1) == as.character(g2)] <- 0)
      
      to_print$time = n
      
      # if in sets there isn't the 0, subset only interesting genes to be printed
      #to_print = to_print[to_print$g2 %in% ge,]
      #to_print = to_print[to_print$g1 %in% gg_in_sets,]
      #to_print = merge(to_print, dl, all = T)
      #to_print = to_print[to_print$type %in% names(genes),]
      #to_print = to_print[to_print$g1 %in% genes,]
      
      to_print$t_hk = ifelse((to_print$g2  %in% true_hk) | (to_print$g1  %in% true_hk), "hk", "n")
      to_print = merge(to_print, dl,by="g2", all.x = T)
      to_print_tot = rbind(to_print_tot, to_print)
    }
  
  
  if (all_conditions) {
    ar.genes = vector()
    for (g_v in unique(to_print_tot$g2[!to_print_tot$g2 %in% unlist(df_genes)])) {
      tmp = to_print_tot[to_print_tot$g2 == g_v,]
      if (nrow(tmp) == length(conditions)*length(genes) & all(tmp$p_val <= p_v)) {
        ar.genes = c(ar.genes,g_v)
      }
    }
    to_print_tot = to_print_tot[to_print_tot$g2 %in% ar.genes,]
  }
  #n = "E60"
  #ldir = "results/2019.08.10/"
  to_print_tot = to_print_tot[order(to_print_tot$p_val, decreasing = F),]  
    

  print(paste("min coex:",min(to_print_tot$coex), "max coex", max(to_print_tot$coex),sep = " "))
  
  #pdf(paste(ldir,"gw_cotan.pdf" ,sep = ""),onefile=T)
  ggplot(to_print_tot, aes(time, factor(g2, 
                                        levels = rev(levels(factor(g2))))))
    geom_tile(aes(fill = coex),colour = "black", show.legend = TRUE) +
    facet_grid( type ~ g1  ,scales = "free", space = "free") + 
    scale_fill_gradient2(low = "#E64B35B2", high = "#8491B4B2", mid = "white", midpoint = 0,na.value = "grey80", space = "Lab", guide = "colourbar", aesthetics = "fill", limits = lim_coex, oob=scales::squish)+ theme(legend.position="bottom")+
    theme(axis.title.x = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.background = element_rect(fill="#8491B44C"),
          strip.text.y  = element_text(size = 9, colour = "#3C5488FF"),
          strip.text.x  = element_text(size = 9,angle= 90 ,colour = "#3C5488FF"),
          axis.title.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.text.y = element_text( size = 9, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),  
          #legend.title = element_blank(),
          legend.text = element_text(color = "#3C5488FF",face ="italic" ),
          legend.position = "bottom",
          legend.title=element_blank(),
          legend.key.height = unit(2, "mm"))#+geom_text(aes(label=ifelse(t_hk == "hk", "H","")), color="grey", size=3)
  
}

#---------------
load_data2 <- function(dir,cond, genes, prefix){
  file_p_val = prefix
  
  #quart = q
  for (n in cond){
    #n = 11
    # get into memory the p-val adjusted for the right quartile
    #assign(paste("p_value_",n,"_adj", sep = ""),read.csv(paste(wdir,dir,"p_value_",n,".csv", sep = ""), header = T, row.names = 1), envir=globalenv())
    
    #assign(paste("p_value_",n,"_adj", sep = ""), fread(paste(wdir,dir,"p_value_",n,".csv", sep = ""), select = 1:1 ))
    #assign(paste("p_value_",n,"_adj", sep = ""), as.data.frame(cbind(get(paste("p_value_",n,"_adj", sep = "")),fread(paste(wdir,dir,"p_value_",n,".csv", sep = ""), select = all_genes ))))
    #assign(paste("p_value_",n,"_adj", sep = ""),column_to_rownames(get(paste("p_value_",n,"_adj", sep = "")),"V1"), envir=globalenv())
    
    #To import the file p_value which means with already apply the Bonferroni correction
    c = Sys.time()
    assign(paste("p_value_",n, sep = ""), fread(paste(dir,file_p_val,n,".csv", sep = ""), select = 1:1 ))
    assign(paste("p_value_",n, sep = ""), as.data.frame(cbind(get(paste("p_value_",n, sep = "")),fread(paste(dir,file_p_val,n,".csv", sep = ""), select = genes ))))
    assign(paste("p_value_",n, sep = ""),column_to_rownames(get(paste("p_value_",n,sep = "")),"V1"), envir=globalenv())
    print(Sys.time()-c)
    # get into memory the cout matrixes for the right quartile
    # this code, using "fread" to read the matrix, allow us to load only the columns containing all the genes that have to be analyzed.
    assign(paste("si_si_mat_",n, sep = ""), fread(paste(dir,"contingency_table_si_si_",n,".csv", sep = ""), select = 1:1 ))
    assign(paste("si_si_mat_",n, sep = ""), as.data.frame(cbind(get(paste("si_si_mat_",n, sep = "")),fread(paste(dir,"contingency_table_si_si_",n,".csv", sep = ""), select = genes ))))
    assign(paste("si_si_mat_",n, sep = ""),column_to_rownames(get(paste("si_si_mat_",n, sep = "")),"V1"), envir=globalenv())
    
    assign(paste("si_no_",n, sep = ""), fread(paste(dir,"contingency_table_si_no_",n,".csv", sep = ""), select = 1:1 ))
    assign(paste("si_no_",n, sep = ""), as.data.frame(cbind(get(paste("si_no_",n, sep = "")),fread(paste(dir,"contingency_table_si_no_",n,".csv", sep = ""), select = genes ))))
    assign(paste("si_no_",n, sep = ""),column_to_rownames(get(paste("si_no_",n, sep = "")),"V1"), envir=globalenv())
    
    assign(paste("no_si_",n, sep = ""), fread(paste(dir,"contingency_table_no_si_",n,".csv", sep = ""), select = 1:1 ))
    assign(paste("no_si_",n, sep = ""), as.data.frame(cbind(get(paste("no_si_",n, sep = "")),fread(paste(dir,"contingency_table_no_si_",n,".csv", sep = ""), select = genes ))))
    assign(paste("no_si_",n, sep = ""),column_to_rownames(get(paste("no_si_",n, sep = "")),"V1"), envir=globalenv())
    
    print(paste("Counts matrix for condition",n,"done. Starting to load the estimated matrixes.",sep = " "))
    
    # get into memory the cout matrixes for the right quartile
    # this code, using "fread" to read the matrix, allow us to load only the columns containing all the genes that have to be analyzed.
    assign(paste("si_si_ex_",n, sep = ""), fread(paste(dir,"estimator_si_si_",n,".csv", sep = ""), select = 1:1 ))
    assign(paste("si_si_ex_",n, sep = ""), as.data.frame(cbind(get(paste("si_si_ex_",n, sep = "")),fread(paste(dir,"estimator_si_si_",n,".csv", sep = ""), select = genes ))))
    assign(paste("si_si_ex_",n, sep = ""),column_to_rownames(get(paste("si_si_ex_",n, sep = "")),"V1"), envir=globalenv())
    
    assign(paste("si_no_ex_",n, sep = ""), fread(paste(dir,"estimator_si_no_",n,".csv", sep = ""), select = 1:1 ))
    assign(paste("si_no_ex_",n, sep = ""), as.data.frame(cbind(get(paste("si_no_ex_",n, sep = "")),fread(paste(dir,"estimator_si_no_",n,".csv", sep = ""), select = genes ))))
    assign(paste("si_no_ex_",n, sep = ""),column_to_rownames(get(paste("si_no_ex_",n, sep = "")),"V1"), envir=globalenv())
    
    assign(paste("no_si_ex_",n, sep = ""), fread(paste(dir,"estimator_no_si_",n,".csv", sep = ""), select = 1:1 ))
    assign(paste("no_si_ex_",n, sep = ""), as.data.frame(cbind(get(paste("no_si_ex_",n, sep = "")),fread(paste(dir,"estimator_no_si_",n,".csv", sep = ""), select = genes ))))
    assign(paste("no_si_ex_",n, sep = ""),column_to_rownames(get(paste("no_si_ex_",n, sep = "")),"V1"), envir=globalenv())
    
    assign(paste("no_no_ex_",n, sep = ""), fread(paste(dir,"estimator_no_no_",n,".csv", sep = ""), select = 1:1 ))
    assign(paste("no_no_ex_",n, sep = ""), as.data.frame(cbind(get(paste("no_no_ex_",n, sep = "")),fread(paste(dir,"estimator_no_no_",n,".csv", sep = ""), select = genes ))))
    assign(paste("no_no_ex_",n, sep = ""),column_to_rownames(get(paste("no_no_ex_",n, sep = "")),"V1"), envir=globalenv())
    
    print(paste("Condition",n,"all done",sep = " "))
  }
}

#-----------------
# it import only p-val and coex matrics for the selected genes (columns). Rows are gw.
load_data3.0 <- function(dir,cond, genes, prefix){
  file_p_val = prefix
  #to add an X in front of the genes starting with a number
  #genes[grep('^([0-9])', genes)] = paste("X",genes[grep('^([0-9])', genes)],sep = "")
  #genes[grep("-",genes)] =  gsub("-",".",genes[grep("-",genes)])
  #quart = q
  for (n in cond){
    
    #To import the file p_value 
    c = Sys.time()
    p_val_temp = fread(paste(dir,file_p_val,n,".csv.gz", sep = ""), select = 1:1 , stringsAsFactors = F)
    p_val_temp = as.data.frame(cbind(p_val_temp,fread(paste(dir,file_p_val,n,".csv.gz", sep = ""), select = genes )))
    p_val_temp = column_to_rownames(p_val_temp,"V1")
    
    #get back to original the column names
    #colnames(p_val_temp)[grep('^X([0-9])', colnames(p_val_temp))] = substring(colnames(p_val_temp)[grep('^X([0-9])', colnames(p_val_temp))], 2)
    
    # to assign the dataframe in the global environment
    assign(paste("p_value_",n, sep = ""),p_val_temp, envir=globalenv())
    
    # coex dataframe
    coex_temp = fread(paste(dir,"coex_",n,".csv.gz", sep = ""), select = 1:1 , stringsAsFactors = F)
    coex_temp = as.data.frame(cbind(coex_temp,fread(paste(dir,"coex_",n,".csv.gz", sep = ""), select = genes )))
    coex_temp = column_to_rownames(coex_temp,"V1")
    
    #get back to original the column names
    #colnames(coex_temp)[grep('^X([0-9])', colnames(coex_temp))] = substring(colnames(coex_temp)[grep('^X([0-9])', colnames(coex_temp))], 2)
    
    # to assign the dataframe in the global environment
    assign(paste("coex_",n, sep = ""),coex_temp, envir=globalenv())
    
    print(Sys.time()-c)
    print(paste("Condition",n," done",sep = " "))
  }
}

#--------------------------
#### Functions to load and extract informations #############
# To load only basic part:  raw and cells
load_basic_files <- function(input_dir,t,s){
  cells = as.data.frame(fread(paste(input_dir,"cells_",t,".csv", sep = "")))
  rownames(cells)=cells$V1
  cells = as.matrix(cells[,2:ncol(cells)])
  real_nu = NA
  print(dim(cells))
  
  print("Cells done!")
  if(t %in% c("11","13","15","17") & s == "M1"){
    print("Loading Mouse 1 data")
    raw = read.csv(paste("../data/2018.12.06/E",t,"5_Only_Cortical_Cells_DGE.txt", sep = ""), sep = "\t", row.names = 1)
    
    colnames(raw)=paste("E",t,"_",colnames(raw), sep = "")
    #raw = raw[rownames(cells),colnames(cells)]
  }else if(t %in% c("E14","P0") & s == "M2"){
    print("Loading mouse 2 data")
    raw = read.csv(paste("../data/2019.08.07/Mouse2/GSE123335_",t,"_combined_matrix.txt.gz",sep = ""),sep = "\t", row.names = 1)
    #raw = raw[rownames(cells),colnames(cells)]
  }else if(t %in% c("E16.5","P0","P05","P18","P23") & s == "H1"){
    print("Loading hippocampus data")
    raw = as.data.frame(fread(paste("../storage/cotan/data/separation_age/",t,"_mouse_dentate_gyrus.csv",sep = "")))
    rownames(raw) = raw$V1
    raw = raw[,2:ncol(raw)]
    colnames(raw)=paste(t,"_",colnames(raw), sep = "")
    #raw = raw[rownames(cells),colnames(cells)]
  }else if(t %in% c("fra1","fra10","fra2") & s == "Sy"){
    print("Loading syntetic data")
    raw = read.csv(paste("../results/2019.10.13/simulated.dataset.",t,"cl.csv",sep = ""), header = T, row.names = 1)
    colnames(raw)=paste(t,"_",colnames(raw), sep = "")
    #raw = raw[rownames(cells),colnames(cells)]
  }else if(t %in% c("sym_E17_cl1_800cs","sym_E17_cl6_800cs","sym_E17_cl6_4000cs","sym_P0_cl15_800cs","sym_P0_cl15_4000cs", "sym_E17_cl1_4000cs","sym_E17_cl1_4000cs") & s == "Sy2"){
    print(paste("Loading syntetic dataset 2.0", t, sep = " "))
    if(t == "sym_E17_cl1_800cs"){
      raw = as.data.frame(fread("../results/2020.02.02/simulated.dataset_E17_sqrt.fra.1cl_800cs.csv.gz",))
      rownames(raw) = raw$V1
      real_nu = raw[1,]
      raw = raw[2:nrow(raw),2:ncol(raw)]
    }else if(t == "sym_E17_cl1_4000cs"){
      raw = as.data.frame(fread("../results/2020.02.02/simulated.dataset_E17_sqrt.fra.1cl_4000cs.csv.gz",))
      rownames(raw) = raw$V1
      real_nu = raw[1,]
      raw = raw[2:nrow(raw),2:ncol(raw)]
    }else if(t == "sym_E17_cl6_800cs"){
      raw = as.data.frame(fread("../results/2020.02.02/simulated.dataset_E17_sqrt.fra.6cl_800cs.csv.gz",))
      rownames(raw) = raw$V1
      real_nu = raw[1,]
      raw = raw[2:nrow(raw),2:ncol(raw)]
    }else if(t == "sym_E17_cl6_4000cs"){
      raw = as.data.frame(fread("../results/2020.02.02/simulated.dataset_E17_sqrt.fra.6cl_4000cs.csv.gz",))
      rownames(raw) = raw$V1
      real_nu = raw[1,]
      raw = raw[2:nrow(raw),2:ncol(raw)]
    }else if(t == "sym_P0_cl15_800cs"){
      raw = as.data.frame(fread("../results/2020.02.02/simulated.dataset_P0_average.fra.15cl_800cs.csv.gz",))
      rownames(raw) = raw$V1
      real_nu = raw[1,]
      raw = raw[2:nrow(raw),2:ncol(raw)]
    }else if(t == "sym_P0_cl15_4000cs"){
      raw = as.data.frame(fread("../results/2020.02.02/simulated.dataset_P0_average.fra.15cl_4000cs.csv.gz",))
      rownames(raw) = raw$V1
      real_nu = raw[1,]
      raw = raw[2:nrow(raw),2:ncol(raw)]
    }
  }else if(t =="CD14" & s == "neg"){
    print("Loading negative dataset")
    raw = Read10X("../data/2020.02.16/filtered_matrices_mex/hg19/")
    pbmc <- CreateSeuratObject(counts = raw, project = t, min.cells = 3, min.features = 50)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc <- subset(pbmc, subset = nFeature_RNA > 50 & nFeature_RNA < 700 & percent.mt < 5 )
    
    
    raw = as.matrix(pbmc[["RNA"]]@counts)
    colnames(raw)=paste(t,"_",colnames(raw), sep = "")
    
  }else if(t =="ercc_f" & s == "neg"){
    print("Loading negative dataset")
    raw = Read10X("../data/2020.02.09/filtered_matrices_mex/ercc92/")
    raw = as.data.frame(raw)
    colnames(raw)=paste(t,"_",colnames(raw), sep = "")
  }
  print(paste("Raw done! Dimension raw: ",dim(raw), sep = ""))
  
  
  return(list(raw, cells, real_nu))
  
}


# To load also the clutered structure
load_files <- function(input_dir,t,s){
  
  #cells = read.csv(paste(input_dir,"cells_",t,".csv", sep = ""), header = T, row.names = 1)
  cells = as.data.frame(fread(paste(input_dir,"cells_",t,".csv", sep = "")))
  rownames(cells)=cells$V1
  cells = as.matrix(cells[,2:ncol(cells)])
  
  print(dim(cells))
  print("Cells done!")
  if(t %in% c("11","13","15","17") & s == "M1"){
    print("Loading Mouse 1 data")
    raw = read.csv(paste("../data/2018.12.06/E",t,"5_Only_Cortical_Cells_DGE.txt", sep = ""), sep = "\t", row.names = 1)
    
    colnames(raw)=paste("E",t,"_",colnames(raw), sep = "")
    #raw = raw[rownames(cells),colnames(cells)]
  }else if(t %in% c("E14","P0") & s == "M2"){
    print("Loading mouse 2 data")
    raw = read.csv(paste("../data/2019.08.07/Mouse2/GSE123335_",t,"_combined_matrix.txt.gz",sep = ""),sep = "\t", row.names = 1)
    #raw = raw[rownames(cells),colnames(cells)]
  }else if(t %in% c("E16.5","P0","P05","P18","P23") & s == "H1"){
    print("Loading hippocampus data")
    raw = read.csv(paste("../rstudio/storage/cotan/data/separation_age/",t,"_mouse_dentate_gyrus.csv",sep = ""), header = T, row.names = 1)
    colnames(raw)=paste(t,"_",colnames(raw), sep = "")
    #raw = raw[rownames(cells),colnames(cells)]
  }else if(t =="CD14" & s == "neg"){
    print("Loading negative dataset")
    raw = Read10X("../data/2020.02.09/filtered_gene_bc_matrices/hg19/")
    raw = as.data.frame(raw)
    colnames(raw)=paste(t,"_",colnames(raw), sep = "")
  }else if(t =="ercc_f" & s == "neg"){
    print("Loading negative dataset")
    raw = Read10X("../data/2020.02.09/filtered_matrices_mex/ercc92/")
    raw = as.data.frame(raw)
    colnames(raw)=paste(t,"_",colnames(raw), sep = "")
  }
  print(paste("Raw done! Dimension raw: ",dim(raw), sep = ""))
  
  #raw <<- raw
  load(paste(input_dir,"PCA_and_tSNE/",t,"_elaborated_data_stucture", sep = ""))
  nu_est = read.csv(paste(input_dir,"nu_est_",t,".csv", sep = ""), row.names = 1)
  print(paste("Nu_est loaded for ",t, "condition. Dimension:",dim(nu_est),sep = ""))
  raw_norm = as.matrix(raw[rownames(cells),colnames(cells)]) %*% diag(1/nu_est[,1]) #[,1]
  colnames(raw_norm)=colnames(cells)
  raw_norm = as.data.frame(raw_norm)
  
  load(file = paste(input_dir,"a_minimization_",t, sep = ""))
  load(paste(input_dir,"mu_estimator",t, sep = ""))
  lambda_i = read.csv(paste(input_dir,"lambda_i_",t,".csv", sep = ""), row.names = 1)
  return(list(raw, cells, nu_est,raw_norm, tot2,mu_estimator,lambda_i, dsq.e14ctx ))
  print(paste("Done sample ",t,sep = ""))
}

# To load only data no cluster structure
load_files2.0 <- function(input_dir,t,s){
  lista =  load_basic_files(input_dir,t,s)
  cells =lista[[2]]
  raw = lista[[1]]
  real_nu = lista[[3]]
  #raw <<- raw
  #load(paste(input_dir,"PCA_and_tSNE/",t,"_elaborated_data_stucture", sep = ""))
  
  nu_est = read.csv(paste(input_dir,"nu_est_",t,".csv", sep = ""), row.names = 1)
  print(paste("Nu_est loaded for ",t, "condition. Dimension:",dim(nu_est),sep = ""))
  raw_norm = as.matrix(raw[rownames(cells),colnames(cells)]) %b/% matrix(nu_est[,1], nrow = 1)
  #colnames(raw_norm)=colnames(cells)
  raw_norm = as.data.frame(raw_norm)
  
  load(file = paste(input_dir,"a_minimization_",t, sep = ""))
  load(paste(input_dir,"mu_estimator",t, sep = ""))
  lambda_i = read.csv(paste(input_dir,"lambda_i_",t,".csv", sep = ""), row.names = 1)
  print(paste("Done sample ",t,sep = ""))
  return(list(raw, cells, nu_est,raw_norm, tot2,mu_estimator,lambda_i,real_nu ))
  
}

count_tables <- function(genes,input_dir,t,s){
  lista = load_basic_files(input_dir,t,s)
  cells =lista[[2]]
  #raw = lista[[1]]
  
  cells = cells[rownames(cells) %in% genes,]
  raw = raw[rownames(raw) %in% genes,]
  
  si_si = as.matrix(cells) %*% t(as.matrix(cells))
  somma = rowSums(cells)
  somma = as.matrix(somma)
  si_any = do.call("cbind", replicate(length(rownames(somma)), somma, simplify = FALSE))
  rm(somma)
  colnames(si_any) = rownames(si_any)
  si_no = si_any - si_si
  si_any = t(si_any)
  no_si = si_any - si_si
  rm(si_any)
  no_no = ncol(cells) - (si_si + no_si + si_no)
  out=list()
  past = vector()
  df2 =  as.data.frame(matrix(ncol = length(genes),nrow =length(genes) ))
  colnames(df2) = genes
  rownames(df2) = genes
  for (g1 in genes) {
    
    for (g2 in genes[!genes %in% c(g1,past)]) {
      #print(paste(g1,g2, sep = ""))
      df = as.data.frame(matrix(nrow = 2,ncol = 2))
      rownames(df) = c(paste(g1,"yes", sep = "_"),paste(g1,"no", sep = "_"))
      colnames(df) = c(paste(g2,"yes", sep = "_"),paste(g2,"no", sep = "_"))
      df[paste(g1,"yes", sep = "_"),paste(g2,"yes", sep = "_")] = si_si[g1,g2]
      df[paste(g1,"yes", sep = "_"),paste(g2,"no", sep = "_")] = si_no[g1,g2]
      df[paste(g1,"no", sep = "_"),paste(g2,"yes", sep = "_")] = no_si[g1,g2]
      df[paste(g1,"no", sep = "_"),paste(g2,"no", sep = "_")] = no_no[g1,g2]
      f =fisher.test(df)
      print(df)
      print(f)
      df2[g1,g2]=f$p.value
      out = c(out,c("data"=list(df)))
    }
    past = c(past,g1)
  }
  lista = list(df2,out)
  return(lista)
}


#old function to compare with Seurat
fisher <- function(g){
  #print(g)
  #print(Sys.time())
  temp_p_val = data.frame(matrix(ncol = 1, nrow = 0))
  #rownames(temp_p_val)= 
  colnames(temp_p_val)=g
  for (j in rownames(si_si)) {
    #   j=  "Satb2"
    col_g_names =c(paste(g,"si", sep = "_"),paste(g,"no", sep = "_") ) 
    row_g_names =c(paste(j,"si",sep = "_"),paste(j,"no",sep = "_"))
    m = as.data.frame(matrix(ncol = 2, nrow = 2))
    colnames(m) = col_g_names
    rownames(m) = row_g_names
    m[paste(j,"si", sep = "_"),paste(g,"si",sep = "_")] = si_si[j,g]
    m[paste(j,"si", sep = "_"),paste(g,"no",sep = "_")] = si_no[j,g]
    m[paste(j,"no", sep = "_"),paste(g,"si",sep = "_")] = no_si[j,g]
    m[paste(j,"no", sep = "_"),paste(g,"no",sep = "_")] = n_cells - sum(m, na.rm = T)
    f =fisher.test(m)
    temp_p_val[j,g] =f$p.value
    #print(temp_p_val[j,g])
  }
  #print(Sys.time() )
  return(temp_p_val)
}


GDI_dataframe <- function(t,dir,s){
  input_dir = dir
  files = load_files2.0(input_dir,t,s)# To load all data regarding the time point
  raw =files[[1]]
  cells=files[[2]] 
  nu_est=files[[3]]
  raw_norm = files[[4]] 
  tot2 =files[[5]]
  mu_estimator = files[[6]] 
  lambda_i = files[[7]]
  #To load the S (statistic) matrix
  S = as.data.frame(fread(paste(input_dir,"S_",t,".csv.gz",sep = "")))
  rownames(S)= S$V1
  S=S[,2:ncol(S)]
  
  
  #GDI = as.array(rowQuantiles(as.matrix(S),probs =0.999 , na.rm = T))
  CD.sorted <- t(apply(t(S),2,sort,decreasing=T))
  CD.sorted = CD.sorted[,1:round(ncol(CD.sorted)/20, digits = 0)]
  CD.sorted = pchisq(as.matrix(CD.sorted), df=1, lower.tail=F)
  
  #CD.sorted <- apply(S,2,sort,decreasing=T)
  #CD.sorted = CD.sorted[1:round(nrow(CD.sorted)/20, digits = 0),]
  #CD.sorted = pchisq(as.matrix(CD.sorted), df=1, lower.tail=F)
  
  GDI = rowMeans(CD.sorted)
  #GDI = colMeans(CD.sorted)
  GDI =as.data.frame(GDI)
  colnames(GDI) = "S"
  
  #GDI$perc.0.1 = as.vector(pchisq(as.matrix(GDI), df=1, lower.tail=F))
  
  sum.raw.norm = log(rowSums(raw_norm))
  
  GDI =  merge(GDI, as.data.frame(sum.raw.norm), by="row.names",all.x=TRUE)
  rownames(GDI)= GDI$Row.names
  
  GDI$log.mean.p = -log(GDI$S)
  GDI$GDI = log(GDI$log.mean.p)
  #quant.p.val[is.infinite(quant.p.val$log.perc.0.1),]$log.perc.0.1 = 1+max(quant.p.val[!is.infinite(quant.p.val$log.perc.0.1),]$log.perc.0.1)
  return(GDI)
}
