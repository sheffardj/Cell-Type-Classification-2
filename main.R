#==============================================================================#
###  Title: AUTO-ANNOTATION OF CELLTYPE: XGBoost & Elastic Net Ensemble
###
###  Final Project: Statistics 454/556 - University of Victoria
### 
###  Contributors: Dayten J. Sheffar  &  Leno S. Rocha  &  Hanan Abousaleh
#
#
#
#
#==============================================================================#
################################### READ ME ####################################
#==============================================================================#
#
#  SUBMITTED: Compute Canada: 24 hours, 12 cpus (1 nodes), 8Gb mem-per-cpu
#  METHOD: 10 Cross Validations of 5-folds. 
#  MODELS: XGBoost with Elastic Net Ensemble
#
#  Estimated Runtime: 20.1 HOURS ------- SUGGESTED RUNTIME: 23:59 HOURS
#  
#  Libraries Used:                                        
#  MASS; class; caret; glmnet; tree; dplyr; ggplot2; RColorBrewer; gridExtra;
#  foreach; doParallel; readr; bigstatsr; rpart; stats; ggpubr; xtable;
#  data.table; reshape2; data.table; xgboost; hrbrthemes;



total_time <- proc.time()[[3]] #time job

set.seed(0) #always :)

print('Load Libraries')
#### Libraries ####
library(MASS)   
library(class)  
library(caret)
library(glmnet)
library(tree)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(foreach)
library(doParallel)
library(readr)
library(bigstatsr)
library(rpart)
library(stats)
library(ggpubr)
library(xtable)
library(data.table)
library(reshape2)
library(data.table)
library(xgboost)
library(hrbrthemes)

print('Load Data')
#### load data ####
dat <- readRDS("final556.rds")
ctype <- as.numeric(as.factor(dat$celltype)); rownames(ctype) <- NULL
condt <- as.numeric(as.factor(dat$condt))-1
count <- t(dat$count); rownames(count) <- NULL
rdf <- cbind(condt, count, ctype) %>% as.data.frame()
rm(dat)
rm(count)

#==============================================================================#
############################## Feature Reduction ###############################
#==============================================================================#

print('Performing Feature Reduction')

# some params
n_classes = length(unique(rdf$ctype))
n_pairs <-  choose(n_classes, 2)
n_comparisons <- (dim(rdf)[2]-1)

# File backed matrices are used in foreach environments!
pvals_tmp  <- as_FBM(matrix(NA, nrow = n_comparisons, ncol = n_pairs))

# indexing function
col_index_fn <- function(class,class2){
  if(class == 0){   ## ONE
    if(class2 == 1){
      return(1)
    }
    if(class2 == 2){
      return(2)
    }
    if(class2 == 3){
      return(3)
    }
    if(class2 == 4){
      return(4)
    }
  }
  if(class == 1){    ##TWO
    if(class2 == 2){
      return(5)
    }
    if(class2 == 3){
      return(6)
    }
    if(class2 == 4){
      return(7)
    }
  }
  if(class == 2){     ##THREE
    if(class2 == 3){
      return(8)
    }
    if(class2 == 4){
      return(9)
    }
  }
  if(class == 3){     ##FOUR
    
    if(class2 == 4){
      return(10)
    }}
}

#### feature reduction finding p-values ####
time_test <- proc.time()[[3]]

# parallelism
print(paste0('On ', 2*detectCores(), ' cores across 2 nodes'))
registerDoParallel(detectCores())
t<-foreach(class = unique(rdf$ctype), .combine='c') %do%{ #i.e., pick "AT1"
  classes <- unique(rdf$ctype)
  remove <- which(classes == class)
  
  foreach(class2 = classes[-c(1:remove)]) %:% #i.e., then do parallely cells 2:5
    foreach( j = 1:n_comparisons) %dopar%{    
      
      tmp1 <- rdf[rdf$ctype == class,]   #subset class 1
      tmp2 <- rdf[rdf$ctype == class2,]  #subset class 2
      tmp1$ctype <- 1 #set class 1 to 1
      tmp2$ctype <- 0 #set class 2 to 0
      tmp <- rbind(tmp1,tmp2) #merge
      
      k <- col_index_fn(class,class2) #column index for FBM
      pvals_tmp[j,k]  <- wilcox.test(tmp[,j]~tmp[,"ctype"])$p.value #mann-whitney test
    }
}
rm(t)
gc()
stopImplicitCluster()

print(paste0('    - ',(proc.time()[[3]]-time_test)/60,' minutes')) #timing

pvals_tmp  <- pvals_tmp[] #overwrite FBM as regular matrix
pvals_tmp[is.na(pvals_tmp)] <- 2 #any errors due to columns of zeroes set >1

print("Selecting Predictors by p-value <= 0.05")
p_min  <- as.matrix(apply(pvals_tmp,  1, min)) %>% as.data.frame() #find row min
p_min$index  <- c(1:(length(p_min[,1]))) #add index so not lose order

p_min2 <- p_min[p_min$V1 != 2,] #remove old NAs overwritten as p-val > 2




png("hist_p_vals.png",width=3000, height = 1200,units='px',res=300)
par(mfrow=c(1,2))
hist(
  p_min[p_min$V1 <= 1.2 , 1],
  breaks = 100,
  xlab = 'p-value',
  main = 'p-values',
  col = 'skyblue3'
)
abline(v = 0.05, col = 'red', lty = 2)
hist(
  p_min[p_min$V1 < 0.05 , 1],
  breaks = 100,
  xlab = 'p-value',
  main = 'p-values <= 0.05',
  col = 'skyblue3')
dev.off()



# there are 896 features below p-value cutoff of 0.05
# condt is the LAST feature included with p-value cutoff 0.05
num_feats_select <- 896

# Remove Zeroes
p_min21 <- p_min[p_min[ ,1]  !=2, ]
p_ord  <- p_min[order(p_min[,1]),  ]
p_best  <- p_ord[1:num_feats_select,] # best num_feats_select without zeroes


df_sml <- rdf[, c(p_best$index,  (n_comparisons+1))] #name dataset to use
fwrite(df_sml, "df_sml.csv")
print('Reduced Dimensions')
print(dim(df_sml))
head(df_sml[,c(1:5,893:897)])

#==============================================================================#
############################### Parameter Init #################################
#==============================================================================#

print('Setting Global Parameters')
n_CVs = 10; nfold = 5
n_models = 1 

num_rows <- dim(df_sml)[1]
num_cols <- dim(df_sml)[2]

idx <- (fread("idx_cv.csv")) %>% as.matrix(); idx <- t(idx)[2:11,]
#idx <- matrix(NA, nrow = n_CVs, ncol = dim(df_sml)[1]) #make n_CVs fold-indexes
#for(ii in 1:n_CVs){idx[ii,] <- createFolds(df_sml$ctype, k = nfold, list=F)} ##ADDED

# required imports to 'foreach'
load_packs <- c('foreach', 'glmnet','MASS','class','rpart','stats', 'xgboost')



#==============================================================================#
############################## Hyperparam Tune #################################
#==============================================================================#

# set xgboost params
xgb_params <- list("objective" = "multi:softmax",
                   "eval_metric" = "mlogloss",
                   "num_class" = 5)

# set grid search params
eta_seq = seq(0.21,0.3,0.005)     # default is .3 # default has been around 0.275
mx_dp_seq = 3                     # default is 6  # discovery suggested using only 3


print("Starting XGBoost Hyperparameter Search")

time_hyper_tune <- proc.time()[[3]]

#### XGBoost ####
# Fit cv.nfold * cv.nround XGB models and save OOF predictions

erro_min=10^5 # initial error set large to be overwritten immediately

# using foreach sequentially
t <- foreach(max_depth = mx_dp_seq, .combine = 'c', .packages = load_packs) %:% 
  foreach(eta = eta_seq, .combine=cbind) %do%{
    
    nround    <- 300 # number of XGBoost rounds to run then run
    cv_model2 <- xgb.cv(params = xgb_params,
                        data=as.matrix(df_sml[,-num_cols]),
                        label= match(df_sml[, num_cols], c(unique(df_sml$ctype))) - 1,
                        nrounds = nround, #early_stopping_rounds = 50,
                        nfold = nfold,
                        eta = eta,
                        max_depth = max_depth,
                        verbose = 0, 
                        prediction = TRUE)
    
    # find test error
    aux = cv_model2$evaluation_log$test_mlogloss_mean
    erro_min_current = min(aux) #then get min test error
    
    #if min error is less than globally set error, write new params to file
    if(erro_min_current < erro_min){
      erro_min <- erro_min_current
      aux2 <- last(cv_model2$evaluation_log$train_mlogloss_mean) #get training loss
      corresp_train_err <- which(aux == erro_min_current)
      c(eta, max_depth, corresp_train_err, aux2)
      
      # eta = learning rate
      # max_depth = number of levels of tree
      # corresp_train_err = nrounds (number trees)
      # aux2 = training error
    }
  }

# timing
print(paste0('Hyperparameter Tuning Time: ', (proc.time()[[3]]-time_hyper_tune)/60,' minutes')) #   
print(paste('t : ', unlist(t)))

entries <- length(unlist(t))
hyperparams <- unlist(t)[(entries - 3):(entries-1)] #retrieve params

#set params
eta <- hyperparams[[1]]
max_depth <- hyperparams[[2]]
nrounds <- hyperparams[[3]]

#write to file
fwrite(list(t),'cv_tuning_df_sml.csv')

#display params to be used in celltype annotation
print(paste0("eta = ", eta,    "
    max_dp = ", max_depth,    "
    nrounds = ", nrounds))

rm(t)
gc()



#==============================================================================#
################################ Model Fitting #################################
#==============================================================================#

# functions to retrieve fold and cross validation indices inside 'foreach'
y.train <- function(fold, cv) {
  match(df_sml[idx[cv, ] != fold, num_cols], c(unique(df_sml$ctype))) - 1
}
y.test  <- function(fold, cv) {
  match(df_sml[idx[cv, ] == fold, num_cols], c(unique(df_sml$ctype))) - 1
}
x.train <- function(fold, cv) {
  as.matrix(df_sml[idx[cv, ] != fold,-num_cols])
}
x.test  <- function(fold, cv) {
  as.matrix(df_sml[idx[cv, ] == fold,-num_cols])
}

# use a filebacked matrix for foreach environment
preds2 <- FBM(ncol = n_models, nrow = n_CVs * dim(df_sml)[1], init = NA)

# set params for model fitting
xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "num_class" = 5, "eta" = eta, "max.depth" = max_depth)

# set the alpha sequence for elastic net for ensemble
alpha_seq <- c(seq(0.1,0.3,by=0.1))


print(paste0('Performing ',nfold,'-fold Cross Validation ', n_CVs, ' times'))

registerDoParallel(detectCores())  # setup parallel backend
time_test <- proc.time()[[3]]

# parallelism
t <- foreach(cv = 1:n_CVs, .combine = 'c', .packages = load_packs) %:%  #NESTED
  foreach(fold = 1:nfold, .combine = 'c') %do% {  #dopar is PARALLEL fold 1,5,2,6,7,3...
    
    # global indexing due to out-of-order predictions due to foreach
    cv_row_set <- 1:num_rows + (cv - 1) * num_rows #i.e. 1:4954 for cv==1
    list_row_set <- cv_row_set[idx[cv,] == fold]   #get fold indices per cv
    
    #============##### XGBoost #####=============#
    # Fit cv.nfold * cv.nround XGB models and save OOF predictions
    xgb_model <- xgboost(params = xgb_params,
                         data = x.train(fold, cv), #check parallel
                         label = y.train(fold, cv), 
                         nrounds = nrounds,
                         verbose = FALSE)
    
    # fits xgboost to get preditions
    y_pred <- matrix(predict(xgb_model, x.test(fold, cv)), ncol=5, byrow=T)
    
    # get celltype from column probabilities
    pr1 <- apply(y_pred, 1, which.max) #predictions
    
    # sort columnwise per row for max
    y_pred_sorted <- t(apply(y_pred, 1, sort, decreasing=T))
    
    # is max col (i.e., col_1) an order of mag bigger than col_2
    id_inc <- which((y_pred_sorted[,2]>round(0.1*y_pred_sorted[,1],1))) 
    
    
    #============##### El.Net #####=============#
    # if more than 10 are poorly predicted then ensemble kicks in
    # with elastic net to fix those poor
    if(length(id_inc)>10){
      time_test <- proc.time()[[3]]
      
      # storage
      ensem_preds <- matrix(NA, ncol=length(alpha_seq), nrow=length(id_inc))
      
      # perform elastic net for relevant poorly predicted xgboost results
      foreach(alph = alpha_seq) %do% { #alpha_seq
        en.fit1 <- cv.glmnet(
          x.train(fold, cv)[], #check parallel = T
          y.train(fold, cv)[],
          alpha = alph,
          family = 'multinomial')
        
        print(paste0("     cv=",cv, "   fold=", fold,"   alpha=",alph)); 
        
        # make new predictions
        en_pred <- as.numeric(predict(en.fit1, x.test(fold, cv), s = 'lambda.1se', type = 'class'))
        ensem_preds[,match(alph, alpha_seq)] <- en_pred[id_inc]
      }
      
      # get class
      p <- apply(ensem_preds, 1, function(x) unique(x)[which.max(tabulate(match(x,unique(x))))])
      pr1[id_inc] <- p+1 
      end_time <- (proc.time()[[3]]-time_test)/60
      print(paste0('    - ',end_time,' minutes')) 
    }
    
    preds2[list_row_set,1] <- pr1  #add to the final predictions matrix (long-ass vector)
  }
rm(t)
end_time <- (proc.time()[[3]]-time_test)/60
print(paste0('Fitting took ',end_time,' minutes'))  # timing
stopImplicitCluster() # close parallel backend


preds2 <- preds2[,1] %>% as.data.frame()
fwrite(preds2,'Ensemb_preds.csv')



#==============================================================================#
################################### SCORING ####################################
#==============================================================================#

print("Beginning Model Scoring")

# storage
auc_models <- matrix(NA, ncol = n_models*5, nrow = n_CVs)
pre_models <- auc_models
rec_models <- auc_models
f11_models <- auc_models

auc_models <- as_FBM(auc_models)
pre_models <- as_FBM(pre_models)
rec_models <- as_FBM(rec_models)
f11_models <- as_FBM(f11_models)


registerDoParallel(detectCores()) # setup parallel backend
tmp <-
  foreach(cv_set = 1:n_CVs, .combine = 'c', .packages = c("foreach")) %:% 
  foreach(model = 1:n_models, .combine = 'c') %dopar% {
    
    cv_row_set <- 1:num_rows + (cv_set - 1) * num_rows   # indexing
    a <- table(df_sml$ctype, preds2[cv_row_set, model])
    
    
    tp.fp <- apply(a,2,sum) ### computing precision, recall, and f1 scores for each cell type
    tp.fn <- apply(a,1,sum) ### https://www.baeldung.com/cs/multi-class-f1-score
    tp <- diag(a)
    prec <- tp/tp.fp
    reca <- tp/tp.fn
    f11 <- 2*prec*reca/(prec+reca) 
    
    
    auc=(1:5)*0 #computing AUC for each cell type
    for (c in 1:5) {
      tmp <- AUC(as.numeric(preds2[cv_row_set,]==c), as.numeric(df_sml$ctype==c)) #bigstatsR
      auc[c] = tmp
    }
    
    auc_models[cv_set, ] <- auc
    pre_models[cv_set, ] <- prec
    rec_models[cv_set, ] <- reca
    f11_models[cv_set, ] <- f11 
    
  }
stopImplicitCluster() # close parallel backend


scores <- cbind( auc_models[], 
                 pre_models[],
                 rec_models[],
                 f11_models[]) %>% as.data.frame()

colnames(scores) <- paste0(rep(c("AUC_", "prec_", "recall_", "f1_"),each=5), 
                           c("B_Cells", "Mesothelial_Cells","Myofibroblasts",
                             "pDCs","Smooth_Muscle_Cells"))
rownames(scores) <- c(paste0('CV', 1:n_CVs))


scores_mean <- cbind("Ens.AUC" = rowMeans(auc_models[]),
                     "Ens.Prec" = rowMeans(pre_models[]),
                     "Ens.Recall" = rowMeans(rec_models[]),
                     "Ens.f1" = rowMeans(f11_models[])) %>% as.data.frame()

rownames(scores_mean) <- c(paste0('CV', 1:n_CVs))

fwrite(scores, paste0('Ensemb_ScoresPerCell.csv')) #save the results
fwrite(scores_mean, paste0('Ensemb_ScoresAveraged.csv'))

print("Scoring Complete")




#==============================================================================#
############################# Signifcance Testing ##############################
#==============================================================================#

#### PERFORMANCE BY CELL TYPE ####
# scores <- fread('Ensemb_ScoresPerCell.csv') %>% as.data.frame()
# scores_mean <- fread('Ensemb_ScoresAveraged.csv') %>% as.data.frame()

# load in benchmark scores as given by Xuekui
scores_eln <- fread('elnet_f1_scores.csv') %>% 
  as.data.frame() %>%
  select(-1)
colnames(scores_eln) <- paste0("eln_",colnames(scores_eln))

# perform wilcox tests PER CELL and stores results
colnameses <- c()
for(col in 1:5){
  assign(paste0('w_',col), wilcox.test(scores[,col+15] - scores_eln[,col], exact = FALSE))
  assign(paste0('bind_',col), cbind(scores[,col+15],scores_eln[,col]))
  colnameses <- c(colnameses, colnames(scores)[col+15], colnames(scores_eln)[col])
}

#### AVERAGED MODEL PERFORMANCE ####

# significance test
wilcox.test(scores_mean[,4] - apply(scores_eln,1,mean),exact=FALSE)


#### PLOTTING ####
# combine results for plotting (col1, col2)
total_mat <- cbind(bind_1, bind_2, bind_3, bind_4, bind_5)
colnames(total_mat) <- colnameses

total_mat <- cbind(total_mat,scores_mean[,4],apply(scores_eln,1,mean))
colnames(total_mat)[11:12] <- c('Total Score','Total Score')
#strip names of '_' and prefixes
nameses <- gsub("f1_","",melt(total_mat)[,2])
nameses <- gsub("eln_","",nameses)
nameses <- gsub("[_]", " ", nameses)
nameses2 <- paste0(nameses,c(rep(", Ensemble",10),rep(", Elastic_Net",10)))

# setup plotting matrix
plot_mat <- as.data.frame(nameses2)
plot_mat$values <-  melt(total_mat)[,3]
plot_mat$group <- nameses
plot_mat$ticks <- rep(c(rep("XGB-EN",10),rep("EN",10)),6)
plot_mat


# ggplot boxplot the scores per cell type
plot_mat %>%
  ggplot(aes(x=as.factor(nameses2), y = values, fill=nameses2)) + 
  geom_boxplot() +
  theme_ipsum() +
  theme(legend.position="none",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.text=element_text(size=9),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c(rep(c("EN","XGB-EN"),5))) +
  facet_wrap(.~group, ncol=6, scales = 'free_x') +
  ylab('F1 Scores') + xlab("") + 
  ggtitle('F1 Scores By Model and Cell Type') +
  scale_fill_brewer(name = 'Legend\nCell_Type.Model', palette = "Paired") -> cell_boxes

cell_boxes

#save the plot to directory
ggsave(
  "box_per_cell.png",
  plot = cell_boxes,
  device = 'png',
  scale = 1,
  width = 14,
  height = 8,
  units = 'in',
  dpi = 300
)



# prepare table for presentation
statistic_value <- c(w_1$statistic, w_2$statistic, w_3$statistic, w_4$statistic, w_5$statistic)
p_values <- c(w_1$p.value, w_2$p.value, w_3$p.value, w_4$p.value, w_5$p.value)
mean_difference <- c(mean(scores[,16] - scores_eln[,1]),
                     mean(scores[,17] - scores_eln[,2]),
                     mean(scores[,18] - scores_eln[,3]),
                     mean(scores[,19] - scores_eln[,4]),
                     mean(scores[,20] - scores_eln[,5]))
mean_values <- c(mean(mean_difference), mean(statistic_value), mean(p_values))
tab1 <- rbind(mean_difference, statistic_value, p_values)
tab2 <- cbind(tab1, mean_values) %>% as.data.frame()
colnames(tab2) <- c(unique(nameses),"averages")
fwrite(tab2, paste0("p_vals_table.csv"))




total_end <- (proc.time()[[3]]-total_time)/60
print(paste0('Total time: ',total_end,' minutes'))  # timing
