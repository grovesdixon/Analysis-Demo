#functions.R


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(ggplot2)
library(DESeq2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(vegan)
library(pheatmap)
library(RColorBrewer)
library(vegan)
library(WGCNA)
library(caret)
library(gtools)
library(glmnet)

# MY FUNCTIONS ------------------------------------------------------------


#function to plot pca from normalized expression counts
build_rld_pca = function(df,
                         coldata,
                         ntop = 25000,
                         pcs = 2){
  #get row varainces and select top
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  #build pca
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d = cbind(data.frame(pca$x[,1:pcs]),
            coldata)
  attr(d, "percentVar") <- percentVar[1:2]
  return(d)
}


#function to plot PCA scatterplot from output from build_rld_pca
plot_rld_pca = function(df,
             group_col = 'treatment',
             pc1 = 1,
             pc2 = 2,
             subtitle = "",
             size = 2,
             legend_title=NULL,
             x_invert=1,
             legend_position = 'none',
             fix_coords = TRUE){
  #select PCs to plot and nvert X if needed
  plt_df = tibble(x = df[,paste('PC', pc1, sep = '')]*x_invert,
                  y = df[,paste('PC', pc2, sep = '')],
                  col = df[,group_col])
  #pull out the percent variances
  percentVar = attr(df, "percentVar")[c(pc1,pc2)]
  #build axis labs
  xlab = paste0(paste0(paste0("PC", pc1), ": "), 
         round(percentVar[pc1] * 100), "%")
  ylab = paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 100), "%")
  g = plt_df %>% 
    ggplot(aes(x=x,
               y=y,
               color=col)) + 
    geom_point(size = size) +
    labs(x = xlab,
         y = ylab,
         subtitle=subtitle,
         color=legend_title) +
    theme(legend.position=legend_position,
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank())
  if (fix_coords){
    g = g + coord_fixed() 
  }
  return(g)
}


#function to plot volcano plot from deseq results
plot_volcano = function(df, TITLE){
  sdf = df %>% 
    arrange(pvalue) %>% 
    mutate(sig = !is.na(padj) & padj < 0.1,
           sig = factor(sig, levels=c(TRUE, FALSE))) %>% 
    as_tibble()
  sdf
  
  g=ggplot(data=sdf,
           aes(x=log2FoldChange,
               y=-log(pvalue, 10),
               color=sig)) +
    geom_point() +
    scale_color_manual(values=c('red', 'black')) + 
    labs(title=TITLE,
         x=bquote(log[2]~difference),
         y=bquote("-"*log[10]~"p"),
         color = 'significant') +
    theme(legend.position='none') +
    lims(x=c(-3,3))
  return(g)
}


#function to read in res files from deseq
read_deseq_res = function(x, n){
  print(paste(n, '...', sep=''))
  load(x)
  d=data.frame(res[,c('log2FoldChange', 'pvalue')])
  colnames(d) = c('lfc', 'p')
  colnames(d) = paste(n, colnames(d), sep='_')
  d$gene = rownames(d)
  return(d)
}



#Function to convert a contingency table from confusion matrix into plotable accuracy percentages
get_acc_pct = function(df, method){
  cm=confusionMatrix(data=df$pred,
                     reference=df$obs,
                     positive='stressed')
  tab=cm$table
  pcts = tab %>% 
    sweep(2,colSums(tab),`/`) %>% 
    data.frame() %>% 
    mutate(Agree=if_else(Prediction==Reference,
                         'Agree',
                         'Disagree'),
           Method = method,
           Reference = if_else(Reference=='stressed',
                               'S',
                               'C'))
  return(list(pcts, cm))
}


#function to load and format prediction stats from random forest and lasso regression
load_pred_stats = function(fileName){
  load(fileName)
  lcmdat = get_acc_pct(lres, 'LR')
  rfcmdat = get_acc_pct(rfres, 'RF')
  lpcts = lcmdat[[1]]
  lcm = lcmdat[[2]]
  rfpcts = rfcmdat[[1]]
  rfcm = rfcmdat[[2]]
  return(list('pcts'=list(lpcts, rfpcts), 'cms'=list(lcm, rfcm)))
}

#build a barplot of prediction stats
plot_pred_stats = function(datList, subtitle){
  dat = datList[["pcts"]] %>% 
    purrr::reduce(rbind)
  cmdat = datList[["cms"]]
  lpval = cmdat[[1]]$overall['AccuracyPValue']
  rfpval = cmdat[[2]]$overall['AccuracyPValue']
  lpsymb = as.character(stars.pval(lpval)[1])
  rfpsymb = as.character(stars.pval(rfpval)[1])
  symdf = tibble(x=c(1,2), y=c(1.01,1.01), lab=c(lpsymb, rfpsymb))
  PSIZE=5
  pnudge=0
  #change capitalization to match figs
  dat$Agree[dat$Agree=='Agree']<-'agree'
  dat$Agree[dat$Agree=='Disagree']<-'disagree'
  dat$Agree = factor(dat$Agree, levels=c('disagree', 'agree'))
  dat %>% 
    group_by(Method, Agree) %>% 
    summarize(pct = sum(Freq)/2) %>% 
    ggplot() +
    geom_bar(aes(x=Method, y=pct, fill=Agree), stat='identity') + 
    geom_text(data=symdf, aes(x=x, y=y, label=lab), size=PSIZE, nudge_y=pnudge) +
    labs(y='Freq', subtitle=subtitle) +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position='none') +
    scale_y_continuous(breaks=c(0,0.5, 1), limits=c(0, 1.1))
}


# COLORS ------------------------------------------------------------------

#get ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols=gg_color_hue(3)
COLOR = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100) # display.brewer.all()
control_col = cols[1]
cluster_a_col = COLOR[70]
cluster_b_col = COLOR[30]
