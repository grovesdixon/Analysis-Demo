#Analysis_Demo.R
rm(list=ls())
source('functions.R')

#This demo concerns the general coral stress response
#The hypothesis is that corals turn their genes up and down
#in a sterotyped manner in response to all forms of stress
#the data are all RNAseq datasets for the genus Acropora on SRA ~6/10/19

# LOOK AT TREATMENT TABLE -------------------------------------------------

#load the trait data from stress studies
ll=load('data/stress_coldata.Rdata')
ll
head(stress_coldata)
dim(stress_coldata)

#build summary table for treatment types
sumTab=stress_coldata %>% 
  group_by(Treatment) %>% 
  summarize(`N projects` = length(unique(my_title)),
            `N samples` = n(),
            `N bleached` = sum(bleached=='yes')) %>% 
  mutate(`N bleached` = ifelse(Treatment=='control',
                               0,
                               `N bleached`))
#check
sumTab


# STRESS PCA ---------------------------------------------------------

#load normalized gene expression controlled for project
ll=load('data/stress_rld.Rdata')
ll
srld_df[1:6,1:6]
dim(srld_df)


#build project PCA
spca_df = build_rld_pca(df = srld_df,
                       coldata = stress_coldata,
                       ntop = 1000,
                       pcs=2)

#plot for stress
stress_plt = plot_rld_pca(spca_df,
                          group_col = 'stress',
                          x_invert=-1,
                          legend_position = 'right')
stress_plt


# PLOT VOLCANO PLOT ----------------------------------------------------------

ll=load('data/stress_deseqResults.Rdata')
ll
head(res)
dim(res)
volc=plot_volcano(data.frame(res), TITLE=NULL) + 
  theme(legend.position = 'right') +
  labs(color='significant')

plot_grid(stress_plt, volc)


# PROJECT HEATMAP ---------------------------------------------------------

#gather the individual paths
deseq_paths = list.files(path = 'deseq_results',
                       pattern = '*deseqResults.Rdata',
                       full.names=TRUE)

#and names
deseq_names0 = list.files(path = 'deseq_results',
                         pattern = '*deseqResults.Rdata',
                         full.names=FALSE)
deseq_names = sub('_deseqResults.Rdata', '', deseq_names0)
names(deseq_paths)=deseq_names
head(deseq_paths)


#read in the files
deseq_list = list()
for (i in 1:length(deseq_paths)){
  x=deseq_paths[i]
  load(x)
  n=deseq_names[i]
  deseq_list[[n]]=read_deseq_res(x,n)
}


#merge into single df
deseq_df = deseq_list %>% 
  purrr::reduce(full_join, by = c('gene'))
rownames(deseq_df)=deseq_df$gene
deseq_df[1:10,1:5]


#isolate log2 fold changes
lfc_df = deseq_df %>% 
  dplyr::select(grep('_lfc', colnames(deseq_df)))
colnames(lfc_df) = sub('_lfc', '', colnames(lfc_df))
lfc_df[1:10,1:5]

#get correlation matrix
c=cor(lfc_df, use="pairwise.complete.obs")
c[1:3,1:3]
dim(c)


#load formatted names for plotting
mod_names = read_excel('data/stress_names_modified.xlsx') %>% 
  column_to_rownames('my_title')
mod_names=mod_names[colnames(c),]
projName = paste(mod_names$Bioproject, mod_names$ref)
treat = mod_names$treat
head(mod_names)


#get hiearchicial cluster
d <- as.dist(1-c)
h=hclust(d, method='average')

#build heatmap
BREAKS=seq(-.1, 0.5, length.out=length(COLOR)+1)
pheatmap(c,
         labels_row = treat,
         labels_col=projName,
         na_col=COLOR[100],
         color=COLOR,
         breaks=BREAKS,
         cluster_rows = h,
         cluster_cols = h,
         fontsize_number = 10,
         treeheight_row=25,
         treeheight_col=0)




# STOPPING POINT ----------------------------------------------------------
#CHECK TIME



# CLASSIFICATION ----------------------------------------------------------

#If there is a general stress response, it should be easy to 
#train a model to identify stressed corals from controls
#did this with lasso logistic regression and random forests


#CLASSIFICATION WITH ALL SAMPLES
ll=load('data/classification_predictions.Rdata')
ll


#look at lasso regression coefficients
cdat$gene=rownames(cdat)
head(cdat)
dim(cdat)
cdat %>% 
  filter(gene !='(Intercept)') %>% 
  mutate(nonZero = coef > 0) %>% 
  pull(nonZero) %>% 
  table()


#look at variable importance
var.imp$gene=rownames(var.imp)
plot(density(var.imp$MeanDecreaseAccuracy))
abline(v=0.001, lty=2, col = 'red')
table(var.imp$MeanDecreaseAccuracy > 0.001)

#look at random forest confusion matrix
confusionMatrix(data=rfres$pred,
                reference=rfres$obs,
                positive="stressed")

#look at category prediction models
pred_res = load_pred_stats('data/classification_predictions.Rdata')
pred_res
pred_bp = plot_pred_stats(pred_res, NULL) + theme(legend.position='top')
pred_bp


# LOOK AT CLUSTERS AMONG RF PREDICTIONS -----------------------------------

#plot tree and choose cutoff
plot(h)
cut.height = 1
abline(h = cut.height, col = "red", lty = 2);

# do assignments
clust = cutreeStatic(h, cutHeight = cut.height, minSize = 4)
table(clust)
clus_a = rownames(c)[(clust==1)]
clus_b = rownames(c)[(clust==2)]

#pull out the samples
clus_a_samples = stress_coldata %>% 
  filter(my_title %in% clus_a) %>% 
  pull(Run)

clus_b_samples = stress_coldata %>% 
  filter(my_title %in% clus_b) %>% 
  pull(Run)

#look at how often the two clusters were mis-identified
head(rfres)
wrong_res = rownames_to_column(rfres, 'Run') %>% 
  mutate(clusA = Run %in% clus_a_samples,
         clusB = Run %in% clus_b_samples,
         wrong = pred!=obs) %>% 
  summarize(A = sum(clusA & wrong)/sum(clusA),
            B = sum(clusB & wrong)/sum(clusB))
wrong_res


#plot barplot of misidentification rate by cluster
wrong_res %>% 
  pivot_longer(A:B,
               names_to = 'cluster',
               values_to = 'misidentification rate') %>% 
  ggplot(aes(x=cluster,
             y=`misidentification rate`)) +
  geom_bar(stat='identity')

