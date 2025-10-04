library(mvtnorm)
library(tidyverse)
library(ggplot2)
library(mombf)
library(stringr)
library(glmnet)
library(ncvreg)
library(parallel)
#BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)


# Helper functions
path = "C:/Users/Usuario/Downloads/PhD-20250908T152636Z-1-001/PhD/Limits on consistent variable selection and external information/Numerical illustrations/"
source(paste0(path,"routines.R"))

# Load data
x = read.table(paste0(path,"colon cancer/tgfb.txt"), header=TRUE)
shortlist= as.character(read.table(paste0(path,"colon cancer/mouse_shortlist.txt"), header=TRUE)[,1])
y = as.matrix(x[,1])
X.design = as.matrix(x[,-1])
inshort= (colnames(X.design) %in% shortlist)
mouselist = which(inshort)
notmouselist = which(!inshort)

# Lists of methods
l0method.vec <- c("kappa.o", "EBIC", "S.EB", "S.A", "S.EB.b",
  "S.A.b")

method.vec <- c(l0method.vec,"lasso.cv", "scad.cv")

set.seed(951)

###################################################

res.sel.l0 <- selectionl0.comp(y,X.design=X.design, block0=mouselist, block1=notmouselist)
res.lasso.scad <- sel.lasso.scad.cv(y,X.design=X.design)
res.sel <- data.frame(res.sel.l0,res.lasso.scad$sel)

res.cvmse.l0 <- cvmse.l0.comp(y, X.design, block0=mouselist, block1=notmouselist, K=10, mc.cores=4)
res.cvmse <- data.frame(res.cvmse.l0$cvmse,res.lasso.scad$cvmse)

result.cvmse.df <- pivot_longer(res.cvmse,cols=method.vec[-1],names_to = 'method',values_to = 'cv.mse')
result.sel.df <- pivot_longer(res.sel[,-1],cols=method.vec[-1],names_to = 'method',values_to = 'sel.model')
result.df <- merge(result.sel.df, result.cvmse.df, by = c('method'))
write.csv(result.df, 'selected.cvmse.csv', row.names = FALSE)


get_annotation_by_block <- function(index.sel.ids,all.ids,mouselist){
  # In mouse list
  index <- index.sel.ids[index.sel.ids%in%mouselist]
  b <- all.ids[index]
  probe_ids= sub("X", "", b)
  annotations= AnnotationDbi::select(hgu133plus2.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID",multiVals='first')
  tab1= cbind(index,annotations)
  rownames(tab1)= NULL
  cat('In mouse list\n')
  print(tab1)
  
  # Not in mouse list 
  index <- index.sel.ids[!index.sel.ids%in%mouselist]
  b <- all.ids[index]
  probe_ids= sub("X", "", b)
  annotations= AnnotationDbi::select(hgu133plus2.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID",multiVals='first')
  annotations <- annotations %>% distinct(PROBEID,.keep_all = TRUE)
  
  tab2= cbind(index,annotations)
  rownames(tab2)= NULL
  cat('Not in mouse list\n') 
  print(tab2)
  res <- rbind(c('Mouse list',NA,NA,NA),tab1,c('Not Mouse list',NA,NA,NA),tab2)
  return(res)
} 

cat('Selected by EBIC\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$EBIC),colnames(X.design),mouselist), 'selected.genes.EBIC.annot.csv', row.names = FALSE)

cat('Selected by S.EB\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$S.EB),colnames(X.design),mouselist), 'selected.genes.S.EB.annot.csv', row.names = FALSE)

cat('Selected by S.EB.b\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$S.EB.b),colnames(X.design),mouselist), 'selected.genes.S.EB.b.annot.csv', row.names = FALSE)

cat('Selected by S.A\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$S.A),colnames(X.design),mouselist), 'selected.genes.S.A.annot.csv', row.names = FALSE)

cat('Selected by S.A.b\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$S.A.b),colnames(X.design),mouselist), 'selected.genes.S.A.b.annot.csv', row.names = FALSE)

cat('Selected by LASSO')
write.csv(get_annotation_by_block(model.char2vec(res.sel$lasso.cv),colnames(X.design),mouselist), 'selected.genes.lasso.annot.csv', row.names = FALSE)

cat('Selected by SCAD')
write.csv(get_annotation_by_block(model.char2vec(res.sel$scad.cv),colnames(X.design),mouselist), 'selected.genes.SCAD.annot.csv', row.names = FALSE)



