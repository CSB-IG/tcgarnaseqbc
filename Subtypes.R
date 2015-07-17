# ssh cristobal@Notron
# ssh -p 5263 cristobal@notron.mine.nu
# cd /home/cristobal/Data/merge
# R

#################################################################################
##Subtypes con pbcmc
#################################################################################
options(width=120)
library("NOISeq")
library("EDASeq")
load(file="FULLGC_FULLLength_TMM.RData")
# load(file="FULLGC_FULLLength_TMM_CPM10_ARSYN.RData")

library("pbcmc")
library("BiocParallel")

##FULLGC_FULLLength_TMM
M<-FULLGC_FULLLength_TMM$M[, FULLGC_FULLLength_TMM$Targets$Group=="E"]
genes<-FULLGC_FULLLength_TMM$Annot[, c("EntrezID", "Symbol.y", "EntrezID")]

##FULLGC_FULLLength_TMM_CPM10_ARSYN
# M<-FULLGC_FULLLength_TMM_CPM10_ARSYN$M[, FULLGC_FULLLength_TMM_CPM10_ARSYN$Targets$Group=="E"]
# genes<-FULLGC_FULLLength_TMM_CPM10_ARSYN$Annot[, c("EntrezID", "Symbol.y", "EntrezID")]

names(genes)<-c("probe", "NCBI.gene.symbol", "EntrezGene.ID")
object<-PAM50(exprs=M, annotation=genes)
object

##Antes ARSyN
# A PAM50 molecular permutation classifier object
# Dimensions:
#             nrow ncol
# exprs      17215  780
# annotation 17215    3
# targets        0    0


##ARSyN
# A PAM50 molecular permutation classifier object
# Dimensions:
#             nrow ncol
# exprs      15281  780
# annotation 15281    3
# targets        0    0

object<-filtrate(object, verbose=TRUE)
object<-classify(object, std="median", verbose=TRUE)
object<-permutate(object, nPerm=10000, pCutoff=0.01, where="fdr",
  corCutoff=0.1, keep=TRUE, seed=1234567890, verbose=TRUE,
  BPPARAM=MulticoreParam(workers=4, progressbar=TRUE))

##FULLGC_FULLLength_TMM
##MEDIAN con cuentas crudas
object
#  Basal   Her2   LumA   LumB Normal 
#    178    113    249    156     84 
table(permutation(object)$subtype$Permuted)
#     Assigned Not Assigned    Ambiguous 
#          435          295           50 

#     Assigned Not Assigned    Ambiguous 
#    55.769231    37.820513     6.410256 

PCA con los Asignados y no Asignados
ARSyN 


Probar que pasa con log(M+1)

object2<-PAM50(exprs=log(M+1), annotation=genes)
object2<-filtrate(object2, verbose=TRUE)
object2<-classify(object2, std="median", verbose=TRUE)
object2<-permutate(object2, nPerm=10000, pCutoff=0.01, where="fdr",
  corCutoff=0.1, keep=TRUE, seed=1234567890, verbose=TRUE,
  BPPARAM=MulticoreParam(workers=4, progressbar=TRUE))

  
object2  
# $subtype
# 
#  Basal   Her2   LumA   LumB Normal 
#    170    115    254    170     71 
table(permutation(object2)$subtype$Permuted)/780*100
# 
#     Assigned    Ambiguous Not Assigned 
#          462           84          234 
# 
#     Assigned    Ambiguous Not Assigned 
#     59.23077     10.76923     30.00000 

##Log(conteos+1)
aux2<-table(permuted=permutation(object2)$subtype$Permuted, PAM50=permutation(object2)$subtype$PAM50)
aux2<-cbind(aux2, Total=rowSums(aux2))
aux2<-rbind(aux2, Total=colSums(aux2))
aux2
#              Basal Her2 LumA LumB Normal Total
# Assigned       142   70  160   59     31   462
# Ambiguous        2   15   32   12     23    84
# Not Assigned    26   30   62   99     17   234
# Total          170  115  254  170     71   780

462/780*100 
# [1] 59.23077
435/780*100
# [1] 55.76923


##Crudos
aux<-table(permuted=permutation(object)$subtype$Permuted, PAM50=permutation(object)$subtype$PAM50)
aux<-cbind(aux, Total=rowSums(aux))
aux<-rbind(aux, Total=colSums(aux))
aux
#              Basal Her2 LumA LumB Normal Total
# Assigned       144   70  141   43     37   435
# Ambiguous        5    6   20    2     17    50
# Not Assigned    29   37   88  111     30   295
# Total          178  113  249  156     84   780

table(permutation(object)$subtype$Permuted)/780*100
#     Assigned Not Assigned    Ambiguous 
#    55.769231    37.820513     6.410256 

##FULLGC_FULLLength_TMM_CPM10_ARSYN.RData  
##NONE
object
# A PAM50 molecular permutation classifier object
# Dimensions:
#            nrow ncol
# exprs        41  780
# annotation   41    3
# targets       0    0
# Classification: 
#             nrow ncol
# probability  780    5
# correlation  780    5
# $subtype
# 
#  Basal   Her2   LumA   LumB Normal 
#      0      2    556    222      0 

# table(permutation(object)$subtype$Permuted)
# 
#     Assigned Not Assigned 
#           35          745 
  
# Obtaining 10000 permutations for 780 subjects... done.
# Error in `rownames<-`(x, value) : 
#   attempt to set 'rownames' on an object with no dimensions
  
##MEDIAN
# $subtype
# 
#  Basal   Her2   LumA   LumB Normal
#    167    144    179    159    131
table(permutation(object)$subtype$Permuted)
# 
#     Assigned Not Assigned    Ambiguous 
#           44          734            2 
##ROBUST
#  Basal   Her2   LumA   LumB Normal 
#    178    130    173    161    138 
#    
table(permutation(object)$subtype$Permuted)
# 
#     Assigned Not Assigned    Ambiguous 
#           46          733            1    





pdf(file="subject1.pdf")  
subjectReport(object, 1)
dev.off()
  




#################################################################################
##3) DIFFERENTIAL EXPRESSION
#################################################################################
