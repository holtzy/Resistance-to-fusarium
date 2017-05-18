	#--------------------------------------------------------------------------------------------------
	#   
	#		ANALYSE des données de FUSARIOSE 2015. Partie Cartographie génétique à partir du mapping sur référence ADr
	#
	#					Script réalisé par 
	#						-Yan Holtz (yan1166@hotmail.com / holtz@supagro.inra.fr)
	#
	#---------------------------------------------------------------------------------------------------



                                                                            
# -------------
# PLAN DU SCRIPT

  1. Préparation des fichiers pour la cartographie (.raw)
    1.1 Préparation pour le merge
    1.2 Merge avec les SNP de la manip capture baites et suppression doublons et allof (R)
    1.3 fin de la préparation pour Carthagène

  2. Première vue sous Carthagène
    2.1 Calculs sous carthagène
    2.2 Repérage des groupes de liaisons (LG)
    2.3 infos sur les groupes de liaisons (R) 

  3. Attribution d un chromosome aux groupes de liaisons                                               
# -------------







#-----------------------------------------------------------------------------------------
########### STEP 1: PREPARATION DES FICHIERS POUR CARTHAGENE ==> (.raw) ##################
#-----------------------------------------------------------------------------------------

# Les reads RNA-SEQ ont été mappés sur la référence ADr. Des SNPs ont été obtenus, puis filtré
# On a obtenu le fichier de SNP suivant:
more /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/SNP/SNP_FUSA_on_ADr_cov10_pval09_clean

# On va travailler dans le dossier carto, ici:
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ_REF_AD



#---------------1.1 Préparation SNP RNA-Seq pour le merge

# Il va falloir merger ces SNPs RNA-Seq avec les SNPs provenant des captures (projet TRAM)
# Pour ce faire plusieurs étapes sont nécessaires.

# Changement format génotype seul :
more /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/SNP/SNP_FUSA_on_ADr_cov10_pval09_clean |  awk '{ toprint = $7"\t"$8 ; for (i = 1; i <= ((NF-8)/2) ; i++){a=9+(i*2)-1 ; a1=substr($a,1,1) ; a2=substr($a,2,1) ; toprint=toprint "\t" a1"/"a2 } print toprint }' | sed 's/>//' | sed 's/-\//-/g' > SNP_FUSA_on_ADr.genot

# I need 122 column, check: 
more SNP_FUSA_on_ADr.genot | awk '{ print NF }' | sort | uniq
# Never more than 64 missing data, check: 
more SNP_FUSA_on_ADr.genot | awk '{ tot=0 ; for (i=3 ; i<=NF ; i=i+1 ) {if($i=="-"){tot+=1} } ; print tot}' | sort | uniq
# Sometimes Dic2 and Silur have the same genotype... 117 times
more SNP_FUSA_on_ADr.genot | awk '{ if($3==$4 && $3!="-"){print $0}}' | wc -l
#Nombre de cas ou les 2 parents sont manquants... 9 times
more SNP_FUSA_on_ADr.genot | awk '{ if($4=="-" && $3=="-"){print $0}}' | wc -l
#Nombre de cas ou un des 2 parent est hétérozygote.. 61 times
more SNP_FUSA_on_ADr.genot | awk '{ a1=substr($3,1,1) ; a2=substr($3,3,1) ; b1=substr($4,1,1) ; b2=substr($4,3,1) ; if((a1!=a2 && a1!="-") || (b1!=b2 && b1!="-")){print $0}}' | wc -l
#Parmi ces cas la, pour combien l'autre parent est missing ? 1 time
more SNP_FUSA_on_ADr.genot | awk '{ a1=substr($3,1,1) ; a2=substr($3,3,1) ; b1=substr($4,1,1) ; b2=substr($4,3,1) ; if((a1!=a2 && a1!="-") || (b1!=b2 && b1!="-")){print $0}}' | awk '{ if($3=="-" || $4=="-")  print $0}' | wc -l

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

# Suppression des individus qui ont trop de données manquantes
a=apply(ALL, 2, function(x){length(which(x=="-"))})
ALL=ALL[ , -which(a/nrow(ALL)*100>90)]

# Suppression d'individus indésirables. (a partir d'une liste, explication dans le fichier excel carto.
SUPPR <- read.table("delete_liste")
for (i in 1:nrow(SUPPR)) {
ALL<-ALL[,!names(ALL)%in%c(as.character(SUPPR[i,]))]
}

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
########### STEP 2: CARTHAGENE: fabrication des LG ##################
#-----------------------------------------------------------------------------------------

module load bioinfo/carthagene/1.2.3
carthagene



# --------- 2.1 Fabrication des LGs (~ 3 minutes)
#-----------------------------------------------------------------------------------------
# A mettre dans un script puis à envoyer en qsub

dsload SNP_BAIT_RNASEQ_SSR_FUSA.raw
mrkmerges
cgout "num_marqueurs.txt"
mrkinfo
cgout ""
cgout "LG.txt"
group 0.2  6.5
cgout ""
cgout "" 
exit

# envoie en qsub
qsub -q normal.q -b yes -cwd -N tmp_LG "module load bioinfo/carthagene/1.2.3 ; carthagene script_LG" 





# --------- 2.1 Récupération des LGs
#-----------------------------------------------------------------------------------------


more num_marqueurs.txt | grep "merges" | grep -v "merged" | sed 's/ [a-z,A-Z].*merges//' | sed 's/  / /g' | sed 's/^ //' > marqueurs_redondants.txt
more num_marqueurs.txt | sed 's/^ //' | sed 's/^ //' | grep "^[0-9]" | sed 's/:.*//' | sed 's/ $//'  | sed 's/ .* / /' > num_marqueurs_clean.txt
more LG.txt | grep ":" | grep "       [0-9]" | sed 's/       //' | awk '{ for (i=3 ; i<=NF ; i++){ print $1"\t"$i} }' > LG_clean
# NOMBRE DE MARQUEURS UNIQUES 
wc -l LG_clean
while read line ; do LG=$(echo $line | cut -f1 -d" ") ; mark=$(echo $line | cut -f2 -d" ") ; echo $line ; all=$(cat marqueurs_redondants.txt | grep "^$mark ") ; for i in $(echo $all) ; do echo $LG" "$i ; done  ; done < LG_clean | sort | uniq > LG_final.txt
rm  LG.txt LG_clean num_marqueurs.txt




# ---------- 2.3 Bilan des LGs (R):
#-----------------------------------------------------------------------------------------

R
LG=read.table("LG_final.txt") ; colnames(LG)=c("LG","mk_ID")

# Barplot du nbr de marqueur par LG
table(LG[,1])

#Nombre de groupe de liaison :
nlevels(as.factor(LG[,1]))
#Nombre de LG avec plus de 30 marqueurs ?
length(table(LG[,1])[table(LG[,1])>30])
#Nombre de LG avec 1 seul marqueur
length(table(LG[,1])[table(LG[,1])==1]) 










#-----------------------------------------------------------------------------------------
########### STEP 3: Attribution d'un chromosome aux groupes de liaisons ##################
#-----------------------------------------------------------------------------------------


# toujours dans le dossier carto
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ_REF_AD


# --------Info sur les LG :
LG=read.table("LG_final.txt") ; colnames(LG)=c("LG","mk_ID")

# Récupération des noms des marqueurs
NUM=read.table("num_marqueurs_clean.txt") ; colnames(NUM)=c("mk_ID" , "mk_name")

# Utilisation de la ref Assaf uniquement pour attribuer chaque LG à un chromosome
POS=read.table("/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/DATA_ASSAF_ISRAEL/hc_genes_info.tab", header=T, sep=",")

# Fichier bilan
data=merge(LG , NUM , by.x=2 , by.y=1 , all.x=T)
#data=data[grep("TRIDC", data$mk_name),]
fun=function(x){strsplit(x,"\\.")[[1]][1] }
data$contig=unlist(lapply(as.character(data[,3]) , fun))
data=merge(data , POS , by.x=4 , by.y=1 , all.x=T)
data=data[order(data[,3]) , ]
data=data[ , c(3:6)]

# Bilan : quel LG va sur quel chromo ?
data$chromo2=substr(data$seqid , 4 , 5)
data$chromo2[ !data$chromo2%in%c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B") ]=NA
bilan=table(data[,1] , factor(data$chromo2) )
nbr_per_LG=table(data[,1])
bilan=cbind( nbr_per_LG , nbr_with_pos=apply(bilan , 1 , sum) , bilan)
bilan=bilan[order(bilan[,2] , decreasing=T), ]

# table d'attribution des chromosomes (tableau vert)
write.table(x=bilan,file="table_attribution_LG_chromosomes", quote=F, col.names=NA , row.names=T, sep="\t")

# On garde les LG avec au moins 5 marqueurs
bilan=bilan[bilan[,2]>=5 , ]  
sum(bilan[,1])

# Voila les groupes a faire (sort un fichier par groupe ! ):
# Dans notre cas c'est très simple, un LG par chromosome.
for(col in c(3:ncol(bilan))){
 file_name=paste("LGs_for_chromosome_",colnames(bilan)[col],sep="")
    list_of_LG=c()
    for(i in c(1:nrow(bilan))){
        if(bilan[i,col]>=3 & bilan[i,col]==max(bilan[i,c(3:ncol(bilan))]) ){
            list_of_LG=c(list_of_LG , rownames(bilan)[i])
            }}
    write.table(file=file_name , list_of_LG , row.names=F , col.names=F , quote=F)
}

#Je range correctement.
mkdir LG_FORMATION
mv LGs* LG_FORMATION/ 












#-----------------------------------------------------------------------------------------
########### STEP 4: Ordonnancement des marqueurs                        ##################
#-----------------------------------------------------------------------------------------

#Tuto carthagene
# sem                        || va prendre l'ordre par défaut des marqueurs et en faire une carte. La carte est ajoutée au "heap"
# nicemapl nicemapd            || méthode euristique, essaie de mettre les forts LOD et les courtes distances ensemble
# mfmapl et mfmapd            || utilise voyaguer du commerce
# heaprint                    || afficher toutes les cartes du heap
# heaprintd                    || afficher toutes les cartes du heap avec beaucoup de détail
# build 10                    || va inverser des marqueurs pour voir si il peut trouver mieux et créer n nouvelles cartes.
# annealing 100 50 0.1 0.8    || essai des trucs pour améliorer...
# greedy 3 1 1 15 0            || idem, essai des trucs.. a fouiller
#Je dois créer des scripts qui appelent du carthagene pour ordonner les groupes de liaison !



##### CARTA_CLASSIQUE

## 
for i in LG_FORMATION/LGs_for_chromosome_* ; do

    chromo=$(echo $i | sed 's/.*_//')
    echo $chromo
    LGs=$(cat $i | tr "\n" " ")
    echo "dsload SNP_BAIT_RNASEQ_SSR_FUSA.raw" > script_ordo_LG_$chromo
    echo "mrkmerges" >> script_ordo_LG_$chromo
    echo "group 0.2  6.5" >> script_ordo_LG_$chromo
    #Je groupe les LG 47 et 36, et ca devient automatiquement la sélection en cours
    echo "groupmerge "$LGs >> script_ordo_LG_$chromo    
    echo "sem" >> script_ordo_LG_$chromo
    echo "nicemapl" >> script_ordo_LG_$chromo
    echo "nicemapd" >> script_ordo_LG_$chromo
    echo "mfmapl" >> script_ordo_LG_$chromo
    echo "mfmapd" >> script_ordo_LG_$chromo
    echo "build 10" >> script_ordo_LG_$chromo
    echo "annealing 100 50 0.1 0.9" >> script_ordo_LG_$chromo
    echo "greedy 10 1 1 15 0" >> script_ordo_LG_$chromo
    echo "cgexport my_map_"$chromo" map "$chromo >> script_ordo_LG_$chromo
    done 
        

chmod 777 script*
for i in $(ls script*); do
qsub -q bigmem.q -b yes -cwd -N cartographie "module load bioinfo/carthagene/1.2.3 ; carthagene $i ;"
done ;




#Transformation du format donné par carthagène
echo -e "group""\t""marker""\t""position" > carte_finale
for i in my_map*xml ; do
chromo=$(cat $i | head -2 | tail -1 | sed 's/.*=\"//' | sed 's/\">//')
cat $i | grep "MARKER" | sed 's/<MARK.*name=\"//' | sed 's/\" position=\"/ /' | sed 's/\">//' | sed 's/^/'$chromo' /' >> carte_finale
##### Il y a un problème avec le xml (manque un marqueur par rapport au fichier classique)
num=$(cat $i | grep "MARKER" | tail -n1 | cut -f 3 -d "=" | sed 's/"//g' | sed 's/>//g')
ecart=$(awk '{print $(NF-1)}' my_map_$chromo)
posi=$(printf %.1f $(echo "$num+$ecart" | bc -l))
marker=$(awk '{print $NF}' my_map_$chromo)
echo -e "$chromo $marker $posi" >> carte_finale
done

rm script_ord* tmp*  my_ma* carto*


# Ajout des marqueurs redondant
while read line ; do mark=$(echo $line | cut -f2 -d" ") ; num=$(cat num_marqueurs_clean.txt | grep $mark$ | cut -f1 -d" ") ; echo $line" "$num ; done < carte_finale  | grep -v "group" > carte_finale2
while read line ; do mark=$(echo $line | cut -f4 -d" ") ; all=$(cat marqueurs_redondants.txt | grep -w $mark)  ; echo $line" "$all ; done < carte_finale2 > carte_finale3
echo -e "group""\t""marker""\t""position" > carte_final.txt
while read line ; do chromo=$(echo $line | cut -f1 -d" ") ; pos=$(echo $line | cut -f3 -d" ") ; all=$(echo $line | cut -f4- -d" ") ; for i in $(echo $all) ; do equ=$(cat num_marqueurs_clean.txt | grep "^$i " | cut -f2 -d" ") ; echo -e $chromo"\t"$equ"\t"$pos  ; done  ;  done < carte_finale3 | sort | uniq >> carte_final.txt
rm num_marqueurs_clean.txt marqueurs_redondants.txt


#RETOURNEMENT DE CHROMOSMES?
for i in 1B 3A 3B 4A 5B 6A 7A ; do
Rscript /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/return_LG_in_map.R carte_final.txt $i ;
done





# ---------------------------------------------------------------------------------------------------------------------------------------------------------










#-----------------------------------------------------------------------------------------
########### STEP 5: AJOUT DES POSITIONS PHYSIQUES PUTATIVES                     ##################
#-----------------------------------------------------------------------------------------
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ_REF_AD

R
library(tidyverse)

# Get the map, and add a column with contig names
map <- read.table("carte_final.txt",header=TRUE)
fun=function(x){strsplit(x,"@")[[1]][1] }
map$contig=unlist(lapply(as.character(map[,2]) , fun))
map$contig=gsub("\\|TRAES.*","",map$contig)
map$contig=gsub("\\..*","",map$contig)

# Correction to the map: I remove 1 marker that is weird..
map=subset(map, marker!="Cluster_10446|Contig2|likelySeq@329")
map$position[ which(map$group=="2A")] = map$position[ which(map$group=="2A")] - 141.3

# Positions physiques de ADr --> simple, il suffit de lire le fichier de la ref dicoccoides associé au fasta
ADr=read.table("/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/DATA_ASSAF_ISRAEL/hc_genes_info.tab", header=T, sep=",")
ADr$seqid=gsub("chr","",ADr$seqid)
ADr=ADr[which(ADr$seqid!="Un"),] %>% droplevels
ADr=ADr %>% mutate(pos=(end+start)/2)
ADr=ADr[ , c(1,2,7)]
colnames(ADr)=c( "contig" ,  "chromo" ,"position" )
ADr$type="ADr"

# Positions physiques de BWr --> il faut blaster BWr sur ADr pour en déduire la position
links=read.table("/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/DATA_ASSAF_ISRAEL/Liaison_BWr_Assaf")
BWr=merge(links, ADr, by.x=2, by.y=1, all=T)
BWr=BWr[ ,2:4]
colnames(BWr)=c( "contig" ,  "chromo" ,"position" )
BWr$type="BWr"

# Positions physiques de EPO --> il faut blaster EPO sur ADr pour en déduite la position!
links=read.table("/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/DATA_ASSAF_ISRAEL/Liaison_EPOr_Assaf")
EPOr=merge(links, ADr, by.x=2, by.y=1, all=T)
EPOr=EPOr[ ,2:4]
colnames(EPOr)=c( "contig" ,  "chromo" ,"position" )
EPOr$type="EPOr"

# Regroupement
POS=rbind(ADr, BWr, EPOr)

# Check si les contig de la map sont dans les contigs des POS -> je dois pas en avoir trop
map$contig[which(!map$contig %in% POS$contig)]

# Merge physique / génétique
data=merge(map, POS, by.x=4, by.y=1 , all.x=T)
data=data[ , -1]
colnames(data) <- c("group","marker","position","group_phy","position_phy","ref")
data=data[order(data$group, data$position) , ]

# Write file
write.table(x=data,file="map_avec_posi_physique.txt",col.names=TRUE,row.names=FALSE,quote=F,sep="\t")

# Faut il retourner des chromosomes?
data %>%
	group_by(group) %>%
	summarize( cor(position, position_phy, use="complete.obs", method="spearman"))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#















#-----------------------------------------------------------------------------------------
########### STEP 6: TRANSFERT DANS DROPBOX                     ##################
#-----------------------------------------------------------------------------------------

# TRANSFERT DANS /DATA
cd ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA 
scp holtz@CC2-login.cirad.fr://gs7k1/projects/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ_REF_AD/map_avec_posi_physique.txt .

# TRANSFERT DANS /SUPPLEMENTARY DATA
cp map_avec_posi_physique.txt ../../SUPPORTING_DATA/

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


