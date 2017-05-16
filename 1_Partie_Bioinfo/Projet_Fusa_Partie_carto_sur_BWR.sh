	#--------------------------------------------------------------------------------------------------
	#   
	#		ANALYSE des données de FUSARIOSE 2015. Partie Cartographie génétique
	#
	#					Script réalisé par 
	#						-Yan Holtz (yan1166@hotmail.com / holtz@supagro.inra.fr)
	#						-Alban Besnard (quasi-uniquement) (albanbesnard@hotmail.fr)
	#
	#---------------------------------------------------------------------------------------------------


#########################################################################################
###                                   SUMMARY                                         ###
#########################################################################################
###################
###################  1. Préparation des fichiers pour la cartographie (.raw)
###################    1.1 Préparation pour le merge
###################    1.2 Merge avec les SNP de la manip capture baites et suppression doublons et allof (R)
###################    1.3 fin de la préparation pour Carthagène
###################
###################  2. Première vue sous Carthagène
###################    2.1 Calculs sous carthagène
###################    2.2 Repérage des groupes de liaisons (LG)
###################    2.3 infos sur les groupes de liaisons (R) 
###################
###################  3. Attribution d'un chromosome aux groupes de liaisons
#########################################################################################







#-----------------------------------------------------------------------------------------
########### STEP 1: PREPARATION DES FICHIERS POUR CARTHAGENE ==> (.raw) ##################
#-----------------------------------------------------------------------------------------

# on part d'un fichier de SNP classique
# il faut aussi le fichier output.alr qui provient de READS2SNP





#---------------1.1 Préparation pour le merge
#-----------------------------------------------------------------------------------------


# Tri sur les individus
head -2 ~/work/7-PROJET_FUSARIOSE/4-READS2SNP/SNP_AVEC_REDONDANCE_ET_ALLOF/output.alr | tail -1 |  sed 's/sp|mapping.//g' | cut -f3- | sed 's/^/SNP\t/' > tmp
cat  ~/work/7-PROJET_FUSARIOSE/4-READS2SNP/SNP_AVEC_REDONDANCE_ET_ALLOF/SNP_super_clean_FUSA_BLE_TENDRE >> tmp
#Changement format génotype seul :
more ~/work/7-PROJET_FUSARIOSE/4-READS2SNP/SNP_AVEC_REDONDANCE_ET_ALLOF/SNP_super_clean_FUSA_BLE_TENDRE |  awk '{ toprint = $7"\t"$8 ; for (i = 1; i <= ((NF-8)/2) ; i++){a=9+(i*2)-1 ; a1=substr($a,1,1) ; a2=substr($a,2,1) ; toprint=toprint "\t" a1"/"a2 } print toprint }' | sed 's/>//' | sed 's/-\//-/g' > SNP_clean_FUSA_BLE_TENDRE.genot

#Vérification : nombre de génotype, nombre de données manquantes max (64 MAX) :
more SNP_clean_FUSA_BLE_TENDRE.genot | awk '{ print NF }' | sort | uniq
more SNP_clean_FUSA_BLE_TENDRE.genot | awk '{ tot=0 ; for (i=3 ; i<=NF ; i=i+1 ) {if($i=="-"){tot+=1} } ; print tot}' | sort | uniq

#Vérification : nbr de cas ou Dic2 et Silur on le meme génotype (et pas missing) #161
more SNP_clean_FUSA_BLE_TENDRE.genot | awk '{ if($3==$4 && $3!="-"){print $0}}' | wc -l
#Nombre de cas ou les 2 parents sont manquants #12
more SNP_clean_FUSA_BLE_TENDRE.genot | awk '{ if($4=="-" && $3=="-"){print $0}}' | wc -l
#Nombre de cas ou un des 2 parent est hétérozygote ? #145
more SNP_clean_FUSA_BLE_TENDRE.genot | awk '{ a1=substr($3,1,1) ; a2=substr($3,3,1) ; b1=substr($4,1,1) ; b2=substr($4,3,1) ; if((a1!=a2 && a1!="-") || (b1!=b2 && b1!="-")){print $0}}' | wc -l
#Parmi ces cas la, pour combien l'autre parent est missing ? #3 pour RNA-seq #7 pour les baites
more SNP_clean_FUSA_BLE_TENDRE.genot | awk '{ a1=substr($3,1,1) ; a2=substr($3,3,1) ; b1=substr($4,1,1) ; b2=substr($4,3,1) ; if((a1!=a2 && a1!="-") || (b1!=b2 && b1!="-")){print $0}}' | awk '{ if($3=="-" || $4=="-")  print $0}' | wc -l

# Transformation format génotypage. A = allèle Dic2, B=Silur. Hétéro = "-". Attention je perds des SNPs qui ont pas les allèles des parents (jen perds 169)
python /NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/programmes/recode_SNP_for_carthagene.py -SNP SNP_clean_FUSA_BLE_TENDRE.genot -out SNP_clean_FUSA_BLE_TENDRE 


# Création d'un fichier avec entete pour faire le merge avec les Baites !
head -2 ~/work/7-PROJET_FUSARIOSE/4-READS2SNP/SNP_AVEC_REDONDANCE_ET_ALLOF/output.alr | tail -1 |  sed 's/sp|mapping.//g' | cut -f5- | sed 's/^/SNP\t/' > SNP_clean_FUSA_tmp 
cat SNP_clean_FUSA_BLE_TENDRE >> SNP_clean_FUSA_tmp 

############################## On va utiliser un fichier similaire des Baites
# Transformation format génotypage. A = allèle Dic2, B=Silur. Hétéro = "-". Attention je perds des SNPs qui ont pas les allèles des parents (jen perds 50)
python /NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/programmes/recode_SNP_for_carthagene.py -SNP /NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/DIC2_SILUR/CAPTURE_ALL_INDIV_07_2014/MAPPING_ON_EPO/SNP/SNP_attendus_AND_bonus_clean.genot -out SNP_clean_BAITES

head -2 /NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/DIC2_SILUR/CAPTURE_ALL_INDIV_07_2014/MAPPING_ON_EPO/READS2SNP/output.alr | tail -1 |  sed 's/sp|resultat_re_mapping.//g' | cut -f5- | sed 's/^/SNP\t/' | sed 's/TT06DC_/DC/g'  > SNP_clean_BAITES_tmp
cat SNP_clean_BAITES >> SNP_clean_BAITES_tmp 

#cleaning
rm tmp

############################## Un dernier pour les SSR
# les SSR sont les 34 derniers du fichiers .raw
# header du nom de "liste_des_genotypes_dans_lordre"

DIRECTORY="/NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/DIC2_SILUR/DART_MICROSAT_AUTRESMARKEURS"
more $DIRECTORY/liste_des_genotypes_dans_le_genotypage_Dart.txt | sed 's/ //g' | sed 's/\./_/g' | /homedir/besnard/prog/awk_transpose_file.sh | sed 's/ /\t/g' | sed 's/^/SNP\t/1' | sed -e 's/\(_[0-9]\t\)/_0\1/g' | sed -e 's/\(_[0-9]$\)/_0\1/g' | sed 's/_0_/_0/g' | sed 's/TT06DC/DC/g' | sed 's/_/./g' > marker_SSR_tmp

tail -n 34 $DIRECTORY/genotypage_Dic2_Silur.raw | sed 's/^*//g' > tmp_marker
# on en garde seulement 21 qui ne crée pas de problème dans ma carte
for i in $(cat liste_SSR_clean) ; do
grep $i tmp_marker >> marker_SSR_tmp ;
done
rm tmp_marker







#-----------------1.2 MERGE SOUS R AVEC SUPPRESSION DES DOUBLONS ET ALLOF
#-----------------------------------------------------------------------------------------


R
# ouverture des fichiers
FUSA<-read.table("SNP_clean_FUSA_tmp",header=TRUE)
BAITES<-read.table("SNP_clean_BAITES_tmp",header=TRUE)
SSR <- read.table("marker_SSR_tmp",header=TRUE)

# merge
ALL<-merge(BAITES,FUSA,all=TRUE)
ALL<-merge(SSR,ALL,all=TRUE)
ALL[is.na(ALL)]<-factor("-")

# Suppression d'individus indésirables. (a partir d'une liste)
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
write.table(x=t(names(ALL)),file="header_SNP_BAITES_FUSA", sep="\t" ,quote=FALSE ,row.names=FALSE ,col.names=FALSE)
quit("no")


######## sous BASH
# On va fusionner les données des SNP_communs aux baites et à la Fusa.
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


#--------------1.3 CREATION du .raw avec les SNP BAITES et FUSA
#-----------------------------------------------------------------------------------------



#Je crée un fichier .raw pour faire une carto avec les SNPs seulement.
echo "data type ri self" > SNP_BAITES_FUSA.raw
a=$(cat merge_SNP_BAITES_FUSA | wc -l)
b=$(cat merge_SNP_BAITES_FUSA | head -1 | awk '{ print NF-1}' )
echo -e $b"\t"$a"\t""0" | sed 's/ //g' | tr "\t" " " >> SNP_BAITES_FUSA.raw
cat  merge_SNP_BAITES_FUSA | tr "\t" " "  | sed 's/Pt\./Pt-/' | sed 's/^/\*/' >>  SNP_BAITES_FUSA.raw
rm merge_SNP_BAITES_FUSA

##################### AJOUT D'individus
# On a ajouté avec les marqueurs de nombreux individus non séquencés ni par les baytes ni par le RNA-seq (32)
# Je les ai laissé pour la suite en espèrant qu'il n'y ai pas de problème (en particulier l'imputation)
# voici la liste:
DC38.14
DC38.22
DC38.31
DC38.40
DC38.48
DC39.12
DC39.22
DC40.09
DC42.20
DC42.27
DC42.29
DC42.30
DC42.31
DC42.32
DC42.33
DC42.34
DC42.37
DC42.42
DC42.47
DC42.50
DC42.51
DC42.52
DC42.54
DC42.56
DC44.12
DC44.18
DC44.23
DC44.29
DC44.04
DC44.40
DC44.44
DC44.09


### cleaning
rm SNP_c* merge* *_tmp



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##











#-----------------------------------------------------------------------------------------
########### STEP 2: Lancement CARTHAGENE ##################
#-----------------------------------------------------------------------------------------

module load bioinfo/carthagene/1.2.3
carthagene

# --------- 2.1 RNA-seq + BAITES (+16000 marqueurs)
#-----------------------------------------------------------------------------------------


dsload /homedir/besnard/work/7-PROJET_FUSARIOSE/B-CARTOGRAPHIE/SNP_BAITES_FUSA.raw
mrkmerges
cgout "num_marqueurs.txt"
mrkinfo
cgout ""
cgout "LG.txt"
group 0.14  7
cgout ""
cgout "" 
exit



# --------- 2.1 RNA-seq + BAITES Récupération des LGs
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


LG=read.table("LG_final.txt") ; colnames(LG)=c("LG","mk_ID")
#Barplot du nbr de marqueur par LG
barplot(table(LG[,1]) , las=2 , ylab="nbr de marqueurs" , xlab="numero du LG") ; abline(h=30)

#Nombre de groupe de liaison :
nlevels(as.factor(LG[,1]))
#Nombre de LG avec plus de 30 marqueurs ?
length(table(LG[,1])[table(LG[,1])>30])
#Nombre de LG avec 1 seul marqueur
length(table(LG[,1])[table(LG[,1])==1]) 






## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##










#-----------------------------------------------------------------------------------------
########### STEP 3: Attribution d'un chromosome aux groupes de liaisons ##################
#-----------------------------------------------------------------------------------------

#### ATTRIUBTION DES LG A DES CHROMOSOMES

#--------Info sur les LG :
PATH=getwd()
LG=read.table("LG_final.txt") ; colnames(LG)=c("LG","mk_ID")
nbr_per_LG=table(LG[,1])

#Récupération des noms des marqueurs
NUM=read.table("num_marqueurs_clean.txt") ; colnames(NUM)=c("mk_ID" , "mk_name")

#Récupération des positions et attributions 
setwd("/NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/BREAD_WHEAT_IWGSC_and_HORDEUM/TRANSCRIT/")
carte_BW=read.table("Physical_map_of_EPO_on_BW.txt" , na.string="-" , header=T , dec=".")
carte_traes=read.table("physical_map_of_BW.txt", na.string="-" , header=T , dec=".")
POS=rbind(carte_BW,carte_traes)

#Fichier bilan
data=merge(LG , NUM , by.x=2 , by.y=1 , all.x=T)
fun=function(x){strsplit(x,"@")[[1]][1] }
data$contig=unlist(lapply(as.character(data[,3]) , fun))
data=merge(data , POS , by.x=4 , by.y=2 , all.x=T)
data=data[order(data[,3]) , ]
data=data[ , c(3:6)]




#Bilan : quel LG va sur quel chromo ?
data$chromo2=substr(data$group , 1 , 2)
data$chromo2[ !data$chromo2%in%c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B") ]=NA
bilan=table(data[,1] , factor(data$chromo2) )
bilan=cbind( nbr_per_LG , nbr_with_pos=apply(bilan , 1 , sum) , bilan)
bilan=bilan[order(bilan[,2] , decreasing=T), ]

setwd(PATH)
### table d'attribution des chromosomes (tableau vert)
write.table(x=bilan,file="table_attribution_LG_chromosomes", quote=F, col.names=NA , row.names=T, sep="\t")


#DIC2 x SILUR
bilan=bilan[bilan[,2]>=5 , ]  
sum(bilan[,1])

#Voila les groupes a faire (sort un fichier par groupe ! ):
#DS
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






## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##










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
    echo "dsload /homedir/besnard/work/7-PROJET_FUSARIOSE/B-CARTOGRAPHIE/SNP_BAITES_FUSA.raw" > script_ordo_LG_$chromo
    echo "mrkmerges" >> script_ordo_LG_$chromo
    echo "group 0.14 7" >> script_ordo_LG_$chromo
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


# DS ONLY --- Ajout des marqueurs redondant !!!!!!
while read line ; do mark=$(echo $line | cut -f2 -d" ") ; num=$(cat num_marqueurs_clean.txt | grep $mark$ | cut -f1 -d" ") ; echo $line" "$num ; done < carte_finale  | grep -v "group" > carte_finale2
while read line ; do mark=$(echo $line | cut -f4 -d" ") ; all=$(cat marqueurs_redondants.txt | grep -w $mark)  ; echo $line" "$all ; done < carte_finale2 > carte_finale3
echo -e "group""\t""marker""\t""position" > carte_final.txt
while read line ; do chromo=$(echo $line | cut -f1 -d" ") ; pos=$(echo $line | cut -f3 -d" ") ; all=$(echo $line | cut -f4- -d" ") ; for i in $(echo $all) ; do equ=$(cat num_marqueurs_clean.txt | grep "^$i " | cut -f2 -d" ") ; echo -e $chromo"\t"$equ"\t"$pos  ; done  ;  done < carte_finale3 | sort | uniq >> carte_final.txt
rm num_marqueurs_clean.txt marqueurs_redondants.txt


#RETOURNEMENT DE CHROMOSMES?
for i in 1B 2A 7B ; do
Rscript /NAS/g2pop/HOLTZ_YAN_DATA/programmes/return_LG_in_map.R carte_final.txt $i ;
done





# --- FABRICATION CARTE AVEC POSITIONS BW :
R
map <- read.table("carte_final.txt",header=TRUE)
DIRECTORY <- getwd()

# la carte en position physique
setwd("/NAS/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/TRANSCRIT/RELEASE_28/")
POS1=read.table("Physical_map_of_EPO_on_BW.txt" , na.string="-" , header=T , dec=".") ; colnames(POS1)=c( "chromo_BW" ,  "contig" ,"position_BW" )
POS2=read.table("physical_map_of_BW.txt", na.string="-" , header=T , dec=".") ;  colnames(POS2)=c( "chromo_BW" ,  "contig" ,"position_BW" )
POS=rbind(POS1,POS2)

#Fichier bilan
POS$chromo_BW=as.factor(substr(POS$chromo_BW , 1 , 2))
fun=function(x){strsplit(x,"@")[[1]][1] }
map$contig=unlist(lapply(as.character(map[,2]) , fun))
data=merge(map , POS , by.x=4 , by.y=2 , all.x=T)
data=data[ , -1]
names(data) <- c("group","marker","posi","group_Americain","posi_physique")
setwd(DIRECTORY)
write.table(x=data,file="map_avec_posi_physique.txt",col.names=TRUE,row.names=FALSE,quote=F,sep="\t")
quit("no") 


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

