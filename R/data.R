##################################
# R source code file for cleaning TCGAOV data
# Git: this is on the TCGAOV repo, master branch
# Created by Sahir,  July 20, 2016
# Updated:
# Notes:
# 
##################################

rm(list = ls())
source("R/packages.R")

# Get Affymetrix Human Genome U133A expression for ovarian cancer patients
u133a.ov <- getTCGA(disease = "OV", data.type = "mRNA_Array", type = "U133", clinical = TRUE)

u133a.ov$clinical %>% dim
u133a.ov$clinical %>% head
u133a.ov$merged.dat %>% dim

u133a.ov$merged.dat[1:10, 1:20]
colnames(u133a.ov$dat)
rownames(u133a.ov$dat)
dim(u133a.ov$dat)

# rows are genes, columns are people
DT <- as.data.table(u133a.ov$dat, keep.rownames = TRUE)
DT[,1,with = F]
setkey(DT, rn)

# downloaded from "http://journals.plos.org/plosone/article/asset?unique&id=info:doi/10.1371/journal.pone.0018064.s015"
# rows are genes
DT.C1 <- as.data.table(read_excel("RawData/Expression/journal.pone.0018064.s015.XLS", sheet = "C1 signature genes"))
setkey(DT.C1, geneSymbol)
DT.C2 <- as.data.table(read_excel("RawData/Expression/journal.pone.0018064.s015.XLS", sheet = "C2 signature genes"))
setkey(DT.C2, geneSymbol)
DT.C4 <- as.data.table(read_excel("RawData/Expression/journal.pone.0018064.s015.XLS", sheet = "C4 signature genes "))
setkey(DT.C4, geneSymbol)
DT.C5 <- as.data.table(read_excel("RawData/Expression/journal.pone.0018064.s015.XLS", sheet = "C5 signature genes"))
setkey(DT.C5, geneSymbol)

DT.C1$geneSymbol %>% unique %>% length()
DT.C2$geneSymbol %>% unique %>% length()
DT.C4$geneSymbol %>% unique %>% length()
DT.C5$geneSymbol %>% unique %>% length()

DT_C1 <- DT[DT.C1]
DT_C1 <- DT_C1[complete.cases(DT_C1)]

DT_C2 <- DT[DT.C1]
DT_C2 <- DT_C1[complete.cases(DT_C1)]

DT_C4 <- DT[DT.C4]
DT_C4 <- DT_C4[complete.cases(DT_C4)]

DT_C5 <- DT[DT.C5]
DT_C5 <- DT_C5[complete.cases(DT_C5)]

indx <- grep('^TCGA', colnames(DT_C1))

for(j in indx){
  set(DT_C1, i=NULL, j=j, value=DT_C1[[j]]*DT_C1[['logFC']])
  set(DT_C2, i=NULL, j=j, value=DT_C2[[j]]*DT_C2[['logFC']])
  set(DT_C4, i=NULL, j=j, value=DT_C4[[j]]*DT_C4[['logFC']])
  set(DT_C5, i=NULL, j=j, value=DT_C5[[j]]*DT_C5[['logFC']])
}

scores <- as.data.table(data.frame(C1 = DT_C1[, colSums(.SD), .SDcols = indx],
C2 = DT_C2[, colSums(.SD), .SDcols = indx],
C4 = DT_C4[, colSums(.SD), .SDcols = indx],
C5 = DT_C5[, colSums(.SD), .SDcols = indx]), keep.rownames = T)


scores[, `:=`(C1 = scale(C1), C2 = scale(C2), C4 = scale(C4), C5 = scale(C5))]
scores[, lapply(.SD, mean), .SDcols = 2:5]
scores[, lapply(.SD, var), .SDcols = 2:5]


scores[, maxCol := which.max(.SD), by = rn, .SDcols = 2:5]
scores[, table(maxCol)]





sum(DT.C1$geneSymbol %in% DT$rn)


glm

?`~`

(form <- a + b ~ c + d)
(LHS <- form[[2]])
(RHS <- form[[3]])

do.call(c, form)

glm

as.formula(LHS)








sum(DT_signature$geneSymbol %in% rownames(u133a.ov$dat))


u133a.ov$merged.dat %>% dim

DT_signature$geneSymbol %>% table
DT_signature$geneSymbol %>% unique %>% length()

DT_signature$ID




DT_pheno <- fread("RawData/Clinical/TCGA Clinical data with TP53 mutation class.csv")
DT_pheno[BCRPATIENTBARCODE %in% "TCGA-61-1895"]


load("RawData/Expression/PT.TCGA_OV.Affy_RMA_DATA.RData")
DT2 <- as.data.table(DATA, keep.rownames = TRUE)
DT2[,1,with = F]

dim(DT)
head(DT)
DT

DT.C1$geneSymbol %>% length()
DT.C1$geneSymbol %>% unique %>% length()

DT_signature$ID %>% unique %>% length()


sum(unique(DT.C1$ID) %in% DT2$rn)

sum(DT$rn %in% DT_signature$ID)


DT3 <- read.table("RawData/Expression/TCGA_489_UE.txt")
head(DT3)
dim(DT3)
rownames(DT3) 

head(DT2)
DT[, 1, with = F]

colnames(DT)

DT_manifest <- fread("RawData/Expression/file_manifest_Expression-Genes.txt")

all(colnames(DT) %in% DT.manifest$`File Name`)


