
	#--------------------------------------------------------------------------------------------------
	#   
	#		ANALYSE des données de FUSARIOSE 2015. Partie BioInformatique.
	#
	#					Script réalisé par 
	#						-Yan Holtz (yan1166@hotmail.com / holtz@supagro.inra.fr)
	#						-Alban Besnard (albanbesnard@hotmail.fr)
	#
	#---------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------  QUE FAIT CE SCRIPT ???

	# La partie bioinfo fait le mapping + détection de SNP
	# La partie carto fait la carte génétique
	# Cette partie crée une matrice de SNP complète, et réalise l'imputation
	
	# En sortie on a une matrice exploitable pour la détection de QTL


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#








#-----------------------------------------------------------------------------------------
########### STEP 1: PREPARATION D'UN FICHIER DE GENOTYPAGE COMPLET ==> (.raw) ##################
#-----------------------------------------------------------------------------------------

# Cette partie ressemble grandement à la partie cartographie. On doit créer un gros fichier de génotypage avec tous les indiv, en RNA seq + en Bait.
# La seule différence avec la partie carto est que cette fois ci on ne supprime aucun individu.

# Les reads RNA-SEQ ont été mappés sur la référence ADr. Des SNPs ont été obtenus, puis filtré
# On a obtenu le fichier de SNP suivant:
more /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/SNP/SNP_FUSA_on_ADr_cov10_pval09_clean

# On va travailler dans le dossier suivant:
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/REF_ADr/SNP



#---------------1.1 Préparation SNP RNA-Seq pour le merge

# Il va falloir merger ces SNPs RNA-Seq avec les SNPs provenant des captures (projet TRAM)
# Pour ce faire plusieurs étapes sont nécessaires.

# Changement format génotype seul :
more /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/SNP/SNP_FUSA_on_ADr_cov10_pval09_clean |  awk '{ toprint = $7"\t"$8 ; for (i = 1; i <= ((NF-8)/2) ; i++){a=9+(i*2)-1 ; a1=substr($a,1,1) ; a2=substr($a,2,1) ; toprint=toprint "\t" a1"/"a2 } print toprint }' | sed 's/>//' | sed 's/-\//-/g' > SNP_FUSA_on_ADr.genot

# Transformation format génotypage. A = allèle Dic2, B=Silur. Hétéro = "-". Attention je perds des SNPs qui ont pas les allèles des parents (jen perds 127)
python /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/recode_SNP_for_carthagene.py -SNP SNP_FUSA_on_ADr.genot -out tmp

# On créé une entête de fichier =  nom des individus dans le bon ordre, et on y ajoute les SNPs
head -2 /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/READS2SNP/output.alr | tail -1 |  sed 's/sp|resultat_mapping.//g' | cut -f5- | sed 's/^/SNP\t/' > SNP_RNASEQ_FUSA_on_ADr 
cat tmp >> SNP_RNASEQ_FUSA_on_ADr 

# Je n'ai plus que 119 colonnes: 120 indiv - 2 pour Dic2 et Silur + 1 pour le nom du SNP
rm tmp *genot




#---------------1.2 Préparation SNP BAITES

# Transformation format génotypage. A = allèle Dic2, B=Silur. Hétéro = "-". Attention je perds des SNPs qui ont pas les allèles des parents (jen perds 55)
python /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/recode_SNP_for_carthagene.py -SNP /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/CAPTURE_ALL_INDIV_07_2014/MAPPING_ON_EPO/SNP/SNP_attendus_AND_bonus_clean.genot -out SNP_clean_BAITES_tmp

# création fichier avec en tete:
head -2 /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/CAPTURE_ALL_INDIV_07_2014/MAPPING_ON_EPO/READS2SNP/output.alr | tail -1 |  sed 's/sp|resultat_re_mapping.//g' | cut -f5- | sed 's/^/SNP\t/' | sed 's/TT06DC_/DC/g'  > SNP_clean_BAITES
cat SNP_clean_BAITES_tmp >> SNP_clean_BAITES

#cleaning
rm *tmp



#---------------1.4 MERGE SOUS R AVEC SUPPRESSION DES DOUBLONS ET ALLOF

R

# ouverture des fichiers
FUSA<-read.table("SNP_RNASEQ_FUSA_on_ADr",header=TRUE)
BAITES<-read.table("SNP_clean_BAITES",header=TRUE)


# merge
ALL<-merge(BAITES,FUSA,all=TRUE)
ALL[is.na(ALL)]<-factor("-")


# REPERAGE des SNPs communs
A<-table(ALL$SNP)
write.table(x=names(A[A==2]),file="SNP_communs.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

## Ecriture de la table
write.table(x=ALL,file="merge_SNP_BAITES_FUSA",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",na="-")

## Ajout d'un header pour la suite
write.table(x=t(names(ALL)),file="header_SNP_BAIT_RNASEQ_SSR_FUSA", sep="\t" ,quote=FALSE ,row.names=FALSE ,col.names=FALSE)
quit("no")

# On va fusionner les données des SNP_communs aux baites et à la Fusa.
# C'est un peu compliqué car si une données est manquante pour une méthode de seq, alors on utilise la donnée de l'autre méthode
cat merge_SNP_BAITES_FUSA > tmp 
LENGTH=$(head -n1 merge_SNP_BAITES_FUSA | wc -w )

for i in $(cat SNP_communs.txt ); do 
a=$(grep $i merge_SNP_BAITES_FUSA)

echo -n $i
for j in $(seq 2 $LENGTH) ; do 
k=$(($j+$LENGTH))
b=$(echo $a | cut -f $j -d" ")
c=$(echo $a | cut -f $k -d" ")
if [ $b = $c ] ; then
echo -ne "\t""$b" ;
elif [ $b = "-" ] ; then
echo  -ne "\t""$c" ;
elif [ $c = "-" ] ; then
echo  -ne "\t""$b" ;
else
echo -ne "\t-" ; # THIS IS NOT SUPPOSED TO HAPPEN (2 fois)
fi ;
done
echo -en "\n"

grep -v $i tmp > tmp2
cat tmp2 > tmp
rm tmp2

done >> tmp3
cat tmp3 >> tmp
rm tmp3
cat tmp > merge_SNP_BAITES_FUSA
rm tmp




#-----------------1.4 CREATION du format .raw

#Je crée un fichier .raw pour faire une carto avec les SNPs seulement.
echo "data type ri self" > SNP_BAIT_RNASEQ_SSR_FUSA.raw
a=$(cat merge_SNP_BAITES_FUSA | wc -l)
b=$(cat merge_SNP_BAITES_FUSA | head -1 | awk '{ print NF-1}' )
echo -e $b"\t"$a"\t""0" | sed 's/ //g' | tr "\t" " " >> SNP_BAIT_RNASEQ_SSR_FUSA.raw
cat  merge_SNP_BAITES_FUSA | tr "\t" " "  | sed 's/Pt\./Pt-/' | sed 's/^/\*/' >>  SNP_BAIT_RNASEQ_SSR_FUSA.raw
rm merge_SNP_BAITES_FUSA

### cleaning
rm SNP_c* merge* *_tmp SNP_RNASEQ_FUSA_on_ADr












#-----------------------------------------------------------------------------------------
########### STEP 2: IMPUTATION ##################
#-----------------------------------------------------------------------------------------


# Copie du fichier .raw# directory
cd  /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/REF_ADr/IMPUTATION/

# Compte de missing data avant imputation? 16239 markers, 43% of missing data
R
data=read.table("/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/SNP/SNP_BAITES_FUSA.raw" , skip=2 , na.strings="-")
my_fun=function(x){ length(x[is.na(x)]) / length(x) * 100  }
a=apply(data[,c(2:ncol(data))],2, my_fun )
#hist(a)
mean(a)



#----------------------------------------------------------------------------------------
######## CREATION DE LA HAPMAP
### Pour cette étape on a besoin d' :
#----1: un fichier .raw (format d'entrée de carthagene)
#----2: un fichier de carte (format chromo'\t'marker'\t'position)
#----3: un header pour le .raw (le header est un fichier txt avec le nom de tous les individus précédés par "SNP" et avec une tabulation entre chaque.)

# rajout du header sur le .raw
cat  ../SNP/header_SNP_BAIT_RNASEQ_SSR_FUSA | sed 's/DC/TT06DC./g'  > tmp_SNP
cat ../SNP/SNP_BAIT_RNASEQ_SSR_FUSA.raw  | sed '1d' | sed '1d' | sed 's/*//g'  >> tmp_SNP

#On merge les marqueurs pour obtenir leurs positions avec les allèles de chaque individus. ! a ne pas perdre de marqueur!!
# 16239 markers dans le fichier génot | 16227 dans la carte
R
A<-read.table("tmp_SNP",header=TRUE)
B<-read.table("/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ_REF_AD/carte_final.txt",header=TRUE)
bilan<-merge(x=B,y=A,by.x=c("marker"),by.y=c("SNP"))

# mise en place du format hapmap
bilan.hap <- data.frame(bilan$marker, "A/B" , bilan$group , bilan$position , NA , NA , NA , NA , NA , NA , NA )
names(bilan.hap) <- c("rs.", "alleles",	"chrom","pos","strand","assembly.","center","protLSID","assayLSID","panelLSID","QCcode")
hapmap<-cbind(bilan.hap,bilan[,-c(1,2,3)])
hapmap$pos<-round(hapmap$pos*10)
write.table(x=hapmap,file="my_hapmap.hmp.txt",quote=F,row.names=F,col.names=TRUE,sep="\t")
q("no")
rm tmp_SNP


# on remplace les tirets et les N des missing values par des N
cat my_hapmap.hmp.txt | sed 's/-/N/g'  > tmp
cat tmp  | sed 's/\<A\>/C/g' | sed 's/\<B\>/G/g' > my_hapmap.hmp.txt
rm tmp


#----------------------------------------------------------------------------------------
######### CREATION DU PEDIGREE
# FAMILLE = "DIC2_BYBLOS" NAME = "DC38.03"  Parent1=dic2
echo -e "Family""\t""Name""\t""Parent1""\t""Parent2""\t""Contribution1""\t""Contribution1""\t""F" > my_pedigree.txt
for i in $(cut -f 2- ../SNP/header_SNP_BAIT_RNASEQ_SSR_FUSA | sed 's/DC/TT06DC./g'); do 
echo -e "DIC2_SILUR\t$i\tdic2\tsilur\t0.5\t0.5\t0.9" >> my_pedigree.txt
done




#----------------------------------------------------------------------------------------
########## PASSAGE DU PLUGIN

module load compiler/gcc/4.9.2
module load bioinfo/vcftools/0.1.12b
module load system/java/jre8
module load bioinfo/tassel/5.2

# SORT de mon hapmap
run_pipeline.pl  -SortGenotypeFilePlugin -inputFile my_hapmap.hmp.txt -outputFile my_hapmap_sorted.hmp.txt -endPlugin

# Passage de l'imputation
run_pipeline.pl -h my_hapmap_sorted.hmp.txt -separate -CallParentAllelesPlugin -p my_pedigree.txt -maxMissing 0.8 -minR 0.2 -logfile My_logfile.txt -bc1 -bcn -w 30 -windowld -endPlugin -ViterbiAlgorithmPlugin -fillgaps true -phet 0.125 -endPlugin -WritePopulationAlignmentPlugin -file imputed_hapmap.hmp.txt -m false -o both -endPlugin
rm *parents.hmp.txt*

zcat imputed_hapmap.hmp.txt.chr1A*nuc* | head -n1  > imputed_hapmap.hmp.txt
for i in $(ls imputed_hapmap.hmp.txt.chr*nuc*); do echo $i ;zcat $i | sed '1d' >> imputed_hapmap.hmp.txt ; rm $i ; done
cat imputed_hapmap.hmp.txt | sed 's/\<N\>/-/g' > tmp
cat tmp  | sed 's/\<C\>/A/g' | sed 's/\<G\>/B/g' > imputed_hapmap.hmp.txt
rm tmp


#----------------------------------------------------------------------------------------
########## FORMAT .csv pour la détection de QTL (sans Dic2 et Silur)
# tous les marqueurs
more imputed_hapmap.hmp.txt|sed 's/\<S\>/-/g' |cut -f 1,12- | /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/awk_transpose_file.sh | sed 's/ /;/g' | sed 's/rs#/SNP/g'   > fichier_genotypage_QTL.csv


#Nbr de missing data restant
R
data=read.table("fichier_genotypage_QTL.csv" , header=T , na.strings="-" , sep=";")
my_fun=function(x){ length(x[is.na(x)]) / length(x) * 100  }
a=apply(data,1, my_fun )
#hist(a)
mean(a)
q("no")
# --> aller vérifier la matrice, l'effet de l'imputation, ce qui se passe avec les indivs pourris plein de données manquantes


#  Copy dans la dropbox fusariose:
cd ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA
scp holtz@CC2-login.cirad.fr://gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/REF_ADr/IMPUTATION/fichier_genotypage_QTL.csv .


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#



























#------------------------------------------------------
# 		STEP 10 : IMPUTATION BWR
#------------------------------------------------------

# Copie du fichier .raw# directory
cd  /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/REF_ADr/IMPUTATION/

# Compte de missing data avant imputation? 16239 markers, 43% of missing data
R
data=read.table("/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/SNP/SNP_BAITES_FUSA.raw" , skip=2 , na.strings="-")
my_fun=function(x){ length(x[is.na(x)]) / length(x) * 100  }
a=apply(data[,c(2:ncol(data))],2, my_fun )
#hist(a)
mean(a)

#----------------------------------------------------------------------------------------
######## CREATION DE LA HAPMAP
### Pour cette étape on a besoin d' :
#----1: un fichier .raw (format d'entrée de carthagene)
#----2: un fichier de carte (format chromo'\t'marker'\t'position)
#----3: un header pour le .raw (le header est un fichier txt avec le nom de tous les individus précédés par "SNP" et avec une tabulation entre chaque.)

# rajout du header sur le .raw
cat  ../SNP/header_SNP_BAITES_FUSA | sed 's/DC/TT06DC./g'  > tmp_SNP
cat ../SNP/SNP_BAITES_FUSA.raw  | sed '1d' | sed '1d' | sed 's/*//g'  >> tmp_SNP

#On merge les marqueurs pour obtenir leurs positions avec les allèles de chaque individus. ! a ne pas perdre de marqueur!!
# 16239 markers dans le fichier génot | 16227 dans la carte
R
A<-read.table("tmp_SNP",header=TRUE)
B<-read.table("/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ/carte_final.txt",header=TRUE)
bilan<-merge(x=B,y=A,by.x=c("marker"),by.y=c("SNP"))

# mise en place du format hapmap
bilan.hap <- data.frame(bilan$marker, "A/B" , bilan$group , bilan$position , NA , NA , NA , NA , NA , NA , NA )
names(bilan.hap) <- c("rs.", "alleles",	"chrom","pos","strand","assembly.","center","protLSID","assayLSID","panelLSID","QCcode")
hapmap<-cbind(bilan.hap,bilan[,-c(1,2,3)])
hapmap$pos<-round(hapmap$pos*10)
write.table(x=hapmap,file="my_hapmap.hmp.txt",quote=F,row.names=F,col.names=TRUE,sep="\t")
q("no")
rm tmp_SNP


# on remplace les tirets et les N des missing values par des N
cat my_hapmap.hmp.txt | sed 's/-/N/g'  > tmp
cat tmp  | sed 's/\<A\>/C/g' | sed 's/\<B\>/G/g' > my_hapmap.hmp.txt
rm tmp


#----------------------------------------------------------------------------------------
######### CREATION DU PEDIGREE
# FAMILLE = "DIC2_BYBLOS" NAME = "DC38.03"  Parent1=dic2
echo -e "Family""\t""Name""\t""Parent1""\t""Parent2""\t""Contribution1""\t""Contribution1""\t""F" > my_pedigree.txt
for i in $(cut -f 2- ../SNP/header_SNP_BAITES_FUSA | sed 's/DC/TT06DC./g'); do 
echo -e "DIC2_SILUR\t$i\tdic2\tsilur\t0.5\t0.5\t0.9" >> my_pedigree.txt
done




#----------------------------------------------------------------------------------------
########## PASSAGE DU PLUGIN

module load compiler/gcc/4.9.2
module load bioinfo/vcftools/0.1.12b
module load system/java/jre8
module load bioinfo/tassel/5.2

# SORT de mon hapmap
run_pipeline.pl  -SortGenotypeFilePlugin -inputFile my_hapmap.hmp.txt -outputFile my_hapmap_sorted.hmp.txt -endPlugin

# Passage de l'imputation
run_pipeline.pl -h my_hapmap_sorted.hmp.txt -separate -CallParentAllelesPlugin -p my_pedigree.txt -maxMissing 0.8 -minR 0.2 -logfile My_logfile.txt -bc1 -bcn -w 30 -windowld -endPlugin -ViterbiAlgorithmPlugin -fillgaps true -phet 0.125 -endPlugin -WritePopulationAlignmentPlugin -file imputed_hapmap.hmp.txt -m false -o both -endPlugin
rm *parents.hmp.txt*

zcat imputed_hapmap.hmp.txt.chr1A*nuc* | head -n1  > imputed_hapmap.hmp.txt
for i in $(ls imputed_hapmap.hmp.txt.chr*nuc*); do echo $i ;zcat $i | sed '1d' >> imputed_hapmap.hmp.txt ; rm $i ; done
cat imputed_hapmap.hmp.txt | sed 's/\<N\>/-/g' > tmp
cat tmp  | sed 's/\<C\>/A/g' | sed 's/\<G\>/B/g' > imputed_hapmap.hmp.txt
rm tmp


#----------------------------------------------------------------------------------------
########## FORMAT .csv pour la détection de QTL (sans Dic2 et Silur)
# tous les marqueurs
more imputed_hapmap.hmp.txt|sed 's/\<S\>/-/g' |cut -f 1,12- | /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/awk_transpose_file.sh | sed 's/ /;/g' | sed 's/rs#/SNP/g'   > fichier_genotypage_QTL.csv


#Nbr de missing data restant
R
data=read.table("fichier_genotypage_QTL.csv" , header=T , na.strings="-" , sep=";")
my_fun=function(x){ length(x[is.na(x)]) / length(x) * 100  }
a=apply(data,1, my_fun )
#hist(a)
mean(a)
q("no")
# --> aller vérifier la matrice, l'effet de l'imputation, ce qui se passe avec les indivs pourris plein de données manquantes

#  Copy dans la dropbox fusariose:
cd ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA
scp holtz@CC2-login.cirad.fr://gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/IMPUTATION/fichier_genotypage_QTL.csv .


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#






