library(data.table)
library(rrBLUP)
library(lme4qtl)
library(lme4)
library(emmeans)
library(foreach)
library(doParallel)
registerDoParallel(50)


# Clean workspace
rm(list=ls())


DIR_TRN="../../Competition_data/Maize_GxE_Competition_Data/Training_Data/"

#------------------------------------
# Training data - pheno
phe_trn=fread(paste0(DIR_TRN,"1_Training_Trait_Data_2014_2021.csv"),
                 data.table=F)
dim(phe_trn) #136012     26
head(phe_trn)
# 1 megagram/hectare = 100 grams/square meter
# Environment  mean  yield  varied  significantly  among  environments, ranging from 5.2 to 13.0 Mg ha^-1
# 5.2 * 100 grams/square meter * 666.7square meters/mu / 1000grams/kg=346.684 kg/mu #1 mu = 666.7 square meters

# K matrix calculated using the dogrm function from: https://github.com/miguelperezenciso/dogrm
# command: grep -v '#' $DIR_TRN/5_Genotype_Data_All_Years_maf_above_pct2.5.recode.vcf | cut -f 10- | sed 's/\// /g'| sed 's/\./9/g' | dogrm -nind 4928 -add -maf 0.05 -maxmiss 0.05 > add_MAF0.025_MISS0.05.G
K = as.matrix(fread('add_MAF0.025_MISS0.05.G',data.table=F))

# genotype data
hybrids = fread('../../Competition_data/Maize_GxE_Competition_Data/Training_Data/out.012.indv',data.table=F,h=F)[,1]
rownames(K) = colnames(K) = hybrids

phe_trn = subset(phe_trn,Hybrid %in% hybrids)

# calculate an inbreeding statistic for each individual as the diagonal of the GRM (ibd from GCTA)
phe_trn$f = diag(K)[match(phe_trn$Hybrid,rownames(K))]-1



Env = unique(phe_trn$Env)[1]
predictions = foreach(Env = unique(phe_trn$Env)) %dopar% {

  print(match(Env,unique(phe_trn$Env)))
  print(paste("#################",Env,"#######################"))
  # ---------------------------------------------------------
  # subset phe_trn for a specific env
  phe_trn_envi=phe_trn[phe_trn$Env==Env,]

  sapply(list(phe_trn_envi$Replicate,
              phe_trn_envi$Block), function(x) unique(x))

  phe_trn_envi$Replicate = factor(phe_trn_envi$Replicate)
  phe_trn_envi$Block = interaction(phe_trn_envi$Replicate,phe_trn_envi$Block)
  if(is.na(phe_trn_envi$Range[1])) phe_trn_envi$Range = 1
  if(is.na(phe_trn_envi$Pass[1])) phe_trn_envi$Pass = 1
  phe_trn_envi$Range = interaction(phe_trn_envi$Replicate,phe_trn_envi$Range)
  phe_trn_envi$Pass = interaction(phe_trn_envi$Replicate,phe_trn_envi$Pass)

  if(length(unique(phe_trn_envi$Replicate))==1) {
    phe_trn_envi$Genotype = phe_trn_envi$Hybrid
    m = relmatLmer(Yield_Mg_ha ~ 1 + f + (1|Block) + (1|Hybrid) + (1|Genotype),phe_trn_envi,relmat = list(Genotype = K))
    gblup = fixef(m)[1] + fixef(m)[2]*diag(K)[match(rownames(ranef(m)$Genotype),rownames(K))]  + as.matrix(t(m@optinfo$relmat$relfac$Genotype) %*% ranef(m)$Genotype[,1])
    U_cond = K[,rownames(gblup)] %*% MASS::ginv(K[rownames(gblup),rownames(gblup)])
    pred = fixef(m)[1] - fixef(m)[2]*diag(K) + U_cond %*% as.matrix(t(m@optinfo$relmat$relfac$Genotype) %*% ranef(m)$Genotype[,1])

  } else {
    phe_trn_envi$Genotype = phe_trn_envi$Hybrid
    m = relmatLmer(Yield_Mg_ha ~ 0+Replicate + f + (1|Block) + (1|Range)+(1|Pass)+(1|Hybrid) + (1|Genotype),phe_trn_envi,relmat = list(Genotype = K))
    fixefs = fixef(m)
    bf = fixefs['f']
    bR = fixefs[grep('Replicate',names(fixefs))]
    gblup = mean(bR) - bf *diag(K)[match(rownames(ranef(m)$Genotype),rownames(K))] + as.matrix(t(m@optinfo$relmat$relfac$Genotype) %*% ranef(m)$Genotype[,1])
    U_cond = K[,rownames(gblup)] %*% MASS::ginv(K[rownames(gblup),rownames(gblup)])
    pred = mean(bR) - bf *diag(K) + U_cond %*% as.matrix(t(m@optinfo$relmat$relfac$Genotype) %*% ranef(m)$Genotype[,1])

  }
  list(
    Env = Env,
    pred = pred
  )
}
names(predictions) = unique(phe_trn$Env)

prediction_matrix = do.call(cbind,lapply(predictions,function(x) x$pred[rownames(K),1]))
colnames(prediction_matrix) = names(predictions)


pred_data = fread('../../Competition_data/Maize_GxE_Competition_Data/Testing_Data/1_Submission_Template_2022.csv',data.table=F)

hybrids_2022 = unique(pred_data$Hybrid)
env_data_full = data.frame(Env = c(colnames(prediction_matrix),unique(pred_data$Env)))
env_data_full$State = substr(env_data_full$Env,1,2)
X_full = model.matrix(~State,env_data_full)
rownames(X_full) = env_data_full$Env
X_train = X_full[colnames(prediction_matrix),]
X_2022 = X_full[unique(pred_data$Env),]
X_pred_2022 = X_2022 %*% MASS::ginv(crossprod(X_train)) %*% t(X_train)

pred_mat_2022 = prediction_matrix[hybrids_2022,rownames(X_train)] %*% t(X_pred_2022)
preds_2022 = reshape2::melt(pred_mat_2022)
rownames(preds_2022) = paste(preds_2022[,1],preds_2022[,2],sep='::')
pred_data$Yield_Mg_ha = preds_2022[paste(pred_data[,2],pred_data[,1],sep='::'),3]

write.csv(pred_data,file = 'DER_direct_GBLUP_final.csv',row.names=F,quote=F)


