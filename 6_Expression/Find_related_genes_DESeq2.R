#Script R

args <- commandArgs(trailingOnly = TRUE)

col_to_use_min=as.numeric(args[1])
col_to_use_max=as.numeric(args[2])


# Threshold of detection:
thres=0.05

#Expression
data <- read.delim("Expression.txt")
rownames(data)=data[,1]
data=data[,-1]
colnames(data)=gsub("DC", "TT06DC." , colnames(data))
rownames(data)=gsub("\\|T.*","",rownames(data) )


#Phenotype matrix:
pheno=read.table("phenotypage_all_fusa.csv" , header=T , sep=";" )
# Keep only the selected phenotypes:
pheno=pheno[ , c(1,c(col_to_use_min:col_to_use_max)[c(col_to_use_min:col_to_use_max) < ncol(pheno)] )]
#numeric
pheno[,-1]=apply(pheno[,-1],2,as.numeric)
# put geno name as rowname
rownames(pheno)=pheno$geno
pheno=pheno[,-1]
# delete columns with only NA
which(apply( pheno , 2 , function(x) all(is.na(x)) )==TRUE)
pheno=pheno[ , ! apply( pheno , 2 , function(x) all(is.na(x)) ) ]



# Library
library(DESeq2)


# A function that compute the DEgenes related to a phenotypic trait
get_DE_genes_from_pheno=function( trait ){

  # TMP On ne prend que les n premières lignes de data
  don=data
  don<-head(data,n=100)

  # sum_expe contains the trait of interest
  sum_expe=data.frame(geno=rownames(pheno),trait=pheno[,trait] )
  sum_expe=na.omit(sum_expe)

  # in the expression matrix, I keep only individuals genotyped for the marker
  don=don[ , which(colnames(don)%in%sum_expe[,1]) ]

  # reorder sum_expe
  sum_expe=sum_expe[match(colnames(don),sum_expe[,1] ), ] 
  rownames(sum_expe)=sum_expe[,1]

  # Call DeSeq2
  dds <- DESeqDataSetFromMatrix(don, sum_expe, formula( ~ trait) )
  dds <- DESeq(dds, test = c("Wald") )
  res <- results(dds)
  return(res)

  # close function
}


# Apply the function to all columns
bilan=data.frame(matrix(0,0,7))
colnames(bilan)=c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","carac")
for(i in 1:ncol(pheno)){
	print(colnames(pheno)[i])
	print(i)
	DE_genes=get_DE_genes_from_pheno(colnames(pheno)[i])
	DE_genes$carac=colnames(pheno)[i]
	res_sig=as.data.frame( DE_genes[ which(DE_genes$padj<thres) , ] )
	bilan=rbind(bilan, res_sig)
}

bilan=data.frame(gene=rownames(bilan), bilan)


# Write the result
name=paste("resultat_DE_pheno_",col_to_use_min,"_to_",col_to_use_max, sep="")
write.table(bilan, file=name, quote=F, row.names=F, col.names=T)











