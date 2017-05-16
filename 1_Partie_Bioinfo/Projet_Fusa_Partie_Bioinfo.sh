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

	# Ce script traite la partie NGS du projet fusariose: récupération des reads - nettoyage - mapping - appel des génotypes - détection et description des SNPs.
	
		#           Les étapes réalisées sont les suivantes : 
		#            0. Raw DATA
		#				0.1 Téléchargements des données sur le site de Toulouse
		#				0.2 Vérification de la qualité des données avec FASTQC
		#
		#            1. Nettoyage des reads
		#            	1.1 CutAdapt
		#            	1.2 filtration sur qualité et taille
		#            	1.3 Appariement des reads (repérage des singletons)
		#				1.4 Création du conf_mapping
		#
		#            2. Mapping sur EPO
		#               2.1 mapping en soi
		#               2.2 Récupération de statistiques
		#               2.3 Parcours des comptages
		#
		#            3. Mapping sur Blé tendre
		#
		#            4. READS2SNP
		#
		#            5. Assemblage Dic2/Silur des reads non mappés
		#
		#            6. Caractérisation des SNPs
		#
		#			 7. Imputation après carto

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#










#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------
# 		STEP 0: RAW DATA
#------------------------------------------------------

# 0/ Localisation du dossier de travail:
# Ce travail a été fait par Alban Besnard lors de son stage de césure.
# Cette partie a donc été réalisée dans son dossier de travail.
cd /homedir/besnard/work/7-PROJET_FUSARIOSE/

# 0.1/ : Je télécharge les données sur le site de toulouse

# Ajout du préfixe du lien
for i in $(cat list_fastq) ; do
echo "http://ng6.toulouse.inra.fr/fileadmin/data/run/5be4abb24/"$i >> list_liens ;
done

# Création d'un dossier raw data
mkdir 0-RAW_DATA
cd 0-RAW_DATA

# Envois des wget depuis les 5 listes différentes (A B C D E)
for i in $(cat LIST_LIENS/list_liens_D) ; do
qsub -q normal.q -b yes -cwd -V -N Download "wget $i" ;
done

# 0.2/ : On vérifie la qualité des données avec un FASTQC
module load system/java/jre8
mkdir FASTQC
for i in ./*R1.fastq* ; do
qsub -q normal.q -V -b yes -cwd -N quality_check "/usr/local/bioinfo/FastQC/default/fastqc -o FASTQC  $i"
done 

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#










#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------
# 		STEP 1: NETTOYAGE DES READS + CONFMAPPING
#------------------------------------------------------


# --- 1.1	Nettoyage par cutadapt 
mkdir CutadaptR1
mkdir CutadaptR2

cd 0-RAW_DATA
#On repère les tags dans les noms de séquences
for i in `ls *R1*` ; do
echo $i
cd /homedir/besnard/work/7-PROJET_FUSARIOSE/0-RAW_DATA
index=${i%_L*}
echo $index
index=${index##*[0-9]_}
echo $index
pat=${i%_R1*}
echo $pat

rev=${i/R1/R2}

#traitement des R1
mkdir ../CutadaptR1/$i
cd ../CutadaptR1/$i
ln -s ../../0-RAW_DATA/$i

# utilisation d'une table pour ajouter les extra-bases aux hexamères
#AGTCAA	CA
#AGTTCC	GT
#ATGTCA	GA
#CCGTCC	CG
#GTCCGC	AC
#GTGAAA	CG
#GTGGCC	TT
#GTTTCG	GA
#CGTACG	TA
#GAGTGG	AT
#ACTGAT	AT
#ATTCCT	TT
bol=$(grep $index ../../table_correspondance | wc -l)
if [ $bol -eq 1 ] ; then
index=$(sed 's/\t//g' ../../table_correspondance | grep $index) 
fi

echo -e ">adapter\nGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"${index}"ATCTCGTATGCCGTCTTCTGCTTG" >adapterP7.fasta 

cd ../..

# traitement des R2
mkdir CutadaptR2/$rev
cd CutadaptR2/$rev
ln -s ../../0-RAW_DATA/$rev
echo -e ">adapter\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" >adapterP5.fasta

echo $index
cd ../..
done

## J'envoie cutadapt_in_chain :
## On supprime les 9 premières bases de chaque séquence (paramètre -u) (imprécision du à la technique de séquençage)
mkdir Cutadapt_out

module load bioinfo/ARCAD/1
qsub -q normal.q -b yes -V -cwd -N Cut_adapt arcad_hts_1_cutadapt_in_chain.pl -i CutadaptR1 -o Cutadapt_out/ -sa adapterP7.fasta -u 9 -noreverse & 
qsub -q normal.q -b yes -V -cwd -N Cut_adapt arcad_hts_1_cutadapt_in_chain.pl -i CutadaptR2 -o Cutadapt_out/ -sa adapterP5.fasta -u 9 -reverse &

# Suppression répertoires Cutadapts
rm -rf CutadaptR1 CutadaptR2





# --- 1.2 Etape de filtration par longueur et qualité
mkdir filtered

cd Cutadapt_out

module load bioinfo/ARCAD/1
for i in `ls *.fastq.gz`; do qsub -q normal.q -b yes -V -cwd -N filter arcad_hts_2_Filter_Fastq_On_Mean_Quality.pl -f=$i -o ../filtered; done

#Suppression répertoire Cutadapt_out
cd ..
rm -r Cutadapt_out





# --- 1.3 repérage des singletons
mkdir 1-CLEAN_DATA

for i in `ls filtered/*R1*`
do pattern=${i%.filtered.fastq.gz};
pattern=${pattern#filtered/};
patternrev=${pattern/R1/R2}
patternsin=${pattern/_R1/};
rev=${i/R1/R2}
module load bioinfo/ARCAD/1
qsub -q normal.q -b yes -V -cwd -N pair perl /usr/local/bioinfo/ARCAD/1/bin/arcad_hts_3_synchronized_paired_fastq.pl -f $i -r $rev -of 1-CLEAN_DATA/${pattern}.cleaned.paired.fastq.gz -or 1-CLEAN_DATA/${patternrev}.cleaned.paired.fastq.gz -os 1-CLEAN_DATA/${patternsin}.single.fastq.gz
done

# information sur les fastq.gz (compte du nbr de reads)
echo -e "individu""\t""nb_reads" > fastq.info
for i in $(ls *fastq.gz) ; do 
a=$(zcat $i | echo $((`wc -l`/4)))
echo -e $i"\t"$a >> fastq.info
done 






# --- 1.4 CREATION du dossier CONF_MAPPING 
mkdir 1-CLEAN_DATA/CONF_MAPPING

# On garde en mémoire l'endroit ou se trouve notre dossier CONF_MAPPING
cd 1-CLEAN_DATA/CONF_MAPPING
chemin=$(pwd) ### 

cd /homedir/besnard/work/7-PROJET_FUSARIOSE/1-CLEAN_DATA/

for i in `ls *R1*paired*` ; do
rev=${i/R1/R2};
here=$(pwd)
sin=$(echo $i | sed 's/_R1.cleaned.paired/.single/')
pat=$(echo $i | sed 's/_R1.*//') ;
tmp=$( echo $pat | cut -f 1 -d "_")

pat=$(grep $tmp $chemin/clefs_individus | cut -f 2)

echo -e "$here/$i\t$pat\tpaired\t$here/$rev">> ${chemin}/conf_mapping;
echo -e "$here/$sin\t$pat\tsingle" >> ${chemin}/conf_mapping;
done

cd ${chemin}

for i in $(cat conf_mapping | awk '{ print $2 }' | sort | uniq) ; do 
cat conf_mapping | grep $i > conf_mapping_$i ;
done

cd ../..

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#









#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------
# 		STEP 2: MAPPING SUR EPO
#------------------------------------------------------

# Partie assez inutile
# En fait on va s'intéresser principalement  au mapping sur le blé tendre qui vient plus tard

# --- 2.1 mapping sur EPO
mkdir 2-MAPPING_SUR_EPO
cd 2-MAPPING_SUR_EPO

### La référence EPO a pour path
ref="/homedir/besnard/RESSOURCES/EPO_106_After_Homeo_Splitter2.fasta" 

#Préparation des scripts de mapping et envoi
for i in $(ls ../1-CLEAN_DATA/CONF_MAPPING/conf_mapping_* | sed 's#../1-CLEAN_DATA/CONF_MAPPING/conf_mapping_##') ; do
echo nohup /usr/local/bioinfo/ARCAD/1/bin/arcad_hts_4_Mapping_Arcad.pl -i ../1-CLEAN_DATA/CONF_MAPPING/conf_mapping_$i -o mapping.$i.bam -r $ref  -mapper bwa -q normal.q \& > script_$i  ; 
done 
chmod 777 script*
for i in sc* ; do ./$i ; done 




# --- 2.2 Récupération d'un bilan mapping et d'un compte de reads

mkdir RESULTS
cd RESULTS
for i in $( ls ../ | grep "bam" | grep -v "metrics" ) ; do
qsub -q normal.q -b yes -cwd -N tmpflagStat "samtools flagstat ../$i > flagstat_$i"
qsub -q normal.q -b yes -cwd -N tmpidxStat "samtools idxstats ../$i > idxstats_$i"
done
#fichier bilan
echo -e "experience""\t""reads_tot""\t""reads_mapped""\t""%" > bilan_mapping ; for i in fl* ; do a=$(echo $i | sed 's/.bam//' | sed 's/.*ing.//')     ; b=$(cat $i | head -1 | awk '{ print $1}' ) ; c=$(cat $i | head -5 | tail -1 | sed 's/ +.*(/\t/' | sed 's/%.*//') ; echo -e $a"\t"$b"\t"$c ; done >>  bilan_mapping
ls idx* > tmp_liste ; /usr/local/bioinfo/python/2.7.9/bin/python2.7 /NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/programmes/Merge_idxstats.py -list tmp_liste

# nettoyage des noms de colonnes et des fichiers temporaires.
sed 's/idxstats_mapping.//g' Bilan_idxstats.txt > tmp.txt
mv tmp.txt Bilan_idxstats.txt
rm flagstat* idx* tmp*





# --- 2.3 Observation des idxstats (sous R)

#IDXstat EPO.
data=read.table("Bilan_idxstats.txt" , header=T)

#Nombre de contig avec au moins 1 read pour chaque indiv ;
#for (i in 2:ncol(data)){a=nrow(data[data[,i]>0 , ]) ; print(colnames(data)[i]) ; print(a) }

#Je récupère les 10 contigs les plus couverts pour chaque expé une par une:
for( i in c(2:3)){
print(colnames(data)[i])
BB_expe=data[,c(1,i)]
BB_expe[,3]=round(BB_expe[,2] / sum(BB_expe[,2]) *100 , 3)
print(head(BB_expe[order(BB_expe[,3] , decreasing=T) , ] , 100))
print(BB_expe[,1])
}

for( i in c(2:121)){
print(colnames(data)[i])
BB_expe=data[,c(1,i)]
BB_expe[,3]=round(BB_expe[,2] / sum(BB_expe[,2]) *100 , 3)
print(head(BB_expe[order(BB_expe[,3] , decreasing=T) , ] , 10))
}

# % de contigs avec au moins x reads au total ?
BB_expe<-data
a=c(1000:300000)
my_fun=function(x){length(BB_expe[BB_expe[,2]>x , 2]) / nrow(BB_expe) *100}
CC=data.frame(nb_reads=a , pourc=unlist(lapply(a ,  my_fun)) )




# --- 2.5 Observation d'un contig en particulier (sous R)

# Contig extrait de l'idxstats (par un grep par exemple)
data=read.table("Contig_aspartic_proteinase")
data=t(data)
plot(data)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#









#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------
# 		STEP 3: MAPPING SUR BLE TENDRE
#------------------------------------------------------

mkdir 3-MAPPING_SUR_BLE_TENDRE
cd 3-MAPPING_SUR_BLE_TENDRE

### ref blé tendre
ref="~/RESSOURCES/transcriptome_BW_without_K_D.fasta"

#Préparation des scripts de mapping et envoi
for i in $(ls ../1-CLEAN_DATA/CONF_MAPPING/conf_mapping_* | sed 's#../1-CLEAN_DATA/CONF_MAPPING/conf_mapping_##') ; do
echo nohup /usr/local/bioinfo/ARCAD/1/bin/arcad_hts_4_Mapping_Arcad.pl -i ../1-CLEAN_DATA/CONF_MAPPING/conf_mapping_$i -o mapping.$i.bam -r $ref  -mapper bwa_mem -q normal.q \& > script_$i  ; 
done 
chmod 777 script*
for i in script_DC40* ; do ./$i ; done 

# Flag and IDXstat
mkdir RESULTS
cd RESULTS
for i in $( ls ../ | grep "bam" | grep -v "metrics" ) ; do
qsub -q normal.q -b yes -cwd -N tmpflagStat "samtools flagstat ../$i > flagstat_$i"
qsub -q normal.q -b yes -cwd -N tmpidxStat "samtools idxstats ../$i > idxstats_$i"
done
#fichier bilan
echo -e "experience""\t""reads_tot""\t""reads_mapped""\t""%" > bilan_mapping ; for i in fl* ; do a=$(echo $i | sed 's/.bam//' | sed 's/.*ing.//')     ; b=$(cat $i | head -1 | awk '{ print $1}' ) ; c=$(cat $i | head -5 | tail -1 | sed 's/ +.*(/\t/' | sed 's/%.*//') ; echo -e $a"\t"$b"\t"$c ; done >>  bilan_mapping
ls idx* > tmp_liste ; /usr/local/bioinfo/python/2.7.9/bin/python2.7 /NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/programmes/Merge_idxstats.py -list tmp_liste

# nettoyage des noms de colonnes et des fichiers temporaires.
sed 's/idxstats_mapping.//g' Bilan_idxstats.txt > tmp.txt
mv tmp.txt Bilan_idxstats.txt
rm flagstat* idx* tmp*

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#








#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------
# 		STEP 3 BIS: MAPPING SUR LES DONNEES DE ASSAF
#------------------------------------------------------

# Cette étape s'effectue dans le dossier:
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD

# La référence diccoides de Assaf utilisée dans cette étude est située ici:
more /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/DATA_ASSAF_ISRAEL/TRIDC_WEWseq_PGSB_20160501_CDS_HighConf_REPR.fasta

# Les reads et le conf mapping sont rangés ici:
ls /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/CLEAN_DATA

#Découpage du confmapping
cp /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/CLEAN_DATA/conf_mapping global_conf_mapping
for i in $(cat global_conf_mapping | awk '{ print $2 }' | sort | uniq) ; do cat global_conf_mapping | grep -P $i"\t" > global_conf_mapping_$i ; done

#Préparation des scripts de mapping et envoi.	
for i in $(ls global_conf_mapping_* | sed 's/global_conf_mapping_//') ; do echo nohup /usr/local/bioinfo/ARCAD/1/bin/arcad_hts_4_Mapping_Arcad.pl -i CONF_MAPPING/global_conf_mapping_$i -o resultat_mapping.$i.bam  -r /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/DATA_ASSAF_ISRAEL/TRIDC_WEWseq_PGSB_20160501_CDS_HighConf_REPR.fasta   -mapper bwa_mem -q normal.q \& > ../script_$i  ; done 
cd .. ; chmod 777 script*

# Envoie sur Dic2 pour commencer tranquillement
./script_silur 
rm script_silur

# Puis envoi sur toute la pop
for i in sc* ; do ./$i ; done

#Je récupère les Flagstats + Idxstats
for i in $( ls ../RESULTAT_MAPPING/ | grep -v "bai" ) ; do 
qsub -q normal.q -b yes -cwd -N tmpflagStat "samtools flagstat ../RESULTAT_MAPPING/$i > flagstat_$i"
qsub -q normal.q -b yes -cwd -N tmpidxStat "samtools idxstats ../RESULTAT_MAPPING/$i > idxstats_$i"
done
#fichier bilan
echo -e "experience""\t""reads_tot""\t""reads_mapped""\t""%" > bilan_mapping ; for i in fl* ; do a=$(echo $i | sed 's/.bam//' | sed 's/.*ing.//')	 ; b=$(cat $i | head -1 | awk '{ print $1}' ) ; c=$(cat $i | head -5 | tail -1 | sed 's/ +.*(/\t/' | sed 's/%.*//') ; echo -e $a"\t"$b"\t"$c ; done >>  bilan_mapping
ls idx* > tmp_liste ; /usr/local/bioinfo/python/2.7.9/bin/python2.7  /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/Merge_idxstats.py -list tmp_liste

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#








#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------
# 		STEP 4: SNP CALLING
#------------------------------------------------------


# Cette étape s'effectue dans le dossier:
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/READS2SNP


#Reads2snp mode classique. Je place Dic2 et Silur devant dans la liste des bams
ls /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/RESULTAT_MAPPING/resultat_mapping.dic2.bam > liste_des_bams
ls /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/RESULTAT_MAPPING/resultat_mapping.silur.bam >> liste_des_bams
ls /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/RESULTAT_MAPPING/resultat_mapping.DC*bam >> liste_des_bams

# Envoi de reads2SNP
qsub -q normal.q -b yes -cwd -N reads2snp "/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/reads2snp.bin  -bamlist liste_des_bams  -bamref   /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/DATA_ASSAF_ISRAEL/TRIDC_WEWseq_PGSB_20160501_CDS_HighConf_REPR.fasta  -fis 0.9 -min 1 -th1 0.1 -par 0 -out output -nbth 40"
#bug - renvoi depuis le ALR
qsub -q bigmem.q -b yes -cwd -N reads2snp "/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/reads2snp.bin  -alr output.alr  -fis 0.9 -min 1 -th1 0.1 -par 0 -out output -nbth 30"


# Puis appel des SNPs
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/SNP
qsub -q normal.q -b yes -cwd -N obtention_snp "/usr/local/bioinfo/python/2.7.9/bin/python2.7 /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/concatenate_alr_and_gen.py  -alr  ../READS2SNP/output.alr  -gen ../READS2SNP/output.gen -out SNP_FUSA_on_ADr_cov10_pval09  -cov 10 -pvalue 0.9"

#Combien de SNP en tout?
wc -l SNP_FUSA_on_ADr_cov10_pval09

# Combien de SNP clean ( FIS>=0.8 / He>=0.32 et >80 indiv) --> 10721
cat SNP_FUSA_on_ADr_cov10_pval09 | awk '{ if ( $5 >= 0.8  && $4>=0.32  && $1>=80 ){print $0}}' | wc -l ;

# Je créé un fichier avec les SNPs clean!
cat SNP_FUSA_on_ADr_cov10_pval09 | awk '{ if ( $5 >= 0.8  && $4>=0.32  && $1>=80 ){print $0}}' > SNP_FUSA_on_ADr_cov10_pval09_clean

--> Utilisation des SNPs pour la carto (voir script carto)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#









#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------
# 		STEP 5: ASSEMBLAGE DES UNMAPPED READS DU MAPPING BWR
#------------------------------------------------------


##### 5.1 récupération des unmapped reads

#Je vais récupérer les unmapped reads bam par bam
for i in $(ls ../3-MAPPING_SUR_BLE_TENDRE/BWA_MEM/*bam) ; do
a=$(echo $i | sed 's/.*ing.//' | sed 's/.bam//')
qsub -q normal.q -b yes -cwd -N get_unmapped_reads "samtools view -f 4 $i | cut -f1 > liste_unmapped_reads_$a"
done

mkdir UNMAPPED_READS
##Puis je vais chercher la séquence et la qualité de ces reads dans le fichier d'origine
for i in  $(cat ../1-CLEAN_DATA/CONF_MAPPING/conf_mapping ) ; do echo $i ; done  | grep "homedir" >liste_fastq
for i in $(tail -n +11 liste_fastq ) ; do
indiv=$(echo $i | sed 's/.*1-CLEAN_DATA\///' | cut -f 1 -d "_")
suffixe=$(echo $i | sed 's/.*1-CLEAN_DATA\///' | cut -f 2 -d "L" | cut -c 4-)
echo $indiv
new_name=$(grep $indiv ~/work/7-PROJET_FUSARIOSE/1-CLEAN_DATA/CONF_MAPPING/clefs_individus | cut -f 2)
qsub -q normal.q -b yes -cwd -N unmapped "/usr/local/bioinfo/python/2.7.9/bin/python2.7 ~/prog/take_some_reads_from_fastq.py -fic $i -liste liste_unmapped_reads_${new_name} -out UNMAPPED_READS/unmapped_reads_${new_name}${suffixe}"
done 

mkdir ASSEMBLAGE_DIC2_SILUR
cd ASSEMBLAGE_DIC2_SILUR
##### 5.2 Assemblage Dic2 Silur
ind=unmapped_reads_silur
qsub -q normal.q -b yes -cwd -V -N Assemblage_abyss "module load bioinfo/abyss/1.5.2 ; abyss-pe name=${ind} k=64 in='../UNMAPPED_READS/${ind}_R1.cleaned.paired.fastq.gz ../UNMAPPED_READS/${ind}_R2.cleaned.paired.fastq.gz' se='../UNMAPPED_READS/${ind}.single.fastq.gz'" ;
ind=unmapped_reads_dic2
qsub -q normal.q -b yes -cwd -V -N Assemblage_abyss "module load bioinfo/abyss/1.5.2 ; abyss-pe name=${ind} k=64 in='../UNMAPPED_READS/${ind}_R1.cleaned.paired.fastq.gz ../UNMAPPED_READS/${ind}_R2.cleaned.paired.fastq.gz' se='../UNMAPPED_READS/${ind}.single.fastq.gz'" ;


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#








#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------
# 		STEP 6: DESCRIPTION DES SNPs OBTENUS sur BWR
#------------------------------------------------------

# Je me place dans le dossier contenant les SNPs:
cd  /NAS/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_BW/SNP 

# Je créé un fichier propre avec l'ordre des individus:
more liste_des_bams | sed 's/.bam//' | sed 's/\/home.*ing.//' > ordre_des_individus

# Passage du script de description des SNPs avant nettoyage
Rscript /NAS/g2pop/HOLTZ_YAN_DATA/programmes/Caracterize_my_SNPs.R   SNP_FUSA_on_BLE_TENDRE_cov10_pval09   description_raw_SNPs.pdf   ordre_des_individus

# Idem pour décrire les SNPs clean!
Rscript /NAS/g2pop/HOLTZ_YAN_DATA/programmes/Caracterize_my_SNPs.R   SNP_super_clean_FUSA_BLE_TENDRE   description_clean_SNPs.pdf   ordre_des_individus

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#









#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------
# 		STEP 6: DESCRIPTION DES SNPs OBTENUS sur ADr
#------------------------------------------------------

# Je me place dans le dossier contenant les SNPs:
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_AD/SNP 

# Je créé un fichier propre avec l'ordre des individus:
more liste_des_bams | sed 's/.bam//' | sed 's/\/gs7.*ing.//' > ordre_des_individus

# Passage du script de description des SNPs avant nettoyage
Rscript /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/Caracterize_my_SNPs.R   SNP_FUSA_on_ADr_cov10_pval09   description_raw_SNPs.pdf   ordre_des_individus

# Idem pour décrire les SNPs clean!
Rscript /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/Caracterize_my_SNPs.R   SNP_FUSA_on_ADr_cov10_pval09_clean   description_clean_SNPs.pdf   ordre_des_individus

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------
# 		STEP 10 : IMPUTATION
#------------------------------------------------------

# directory
cd  /NAS/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/IMPUTATION

# Compte de missing data avant imputation? 16239 markers, 43% of missing data
R
data=read.table("/NAS/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/SNP/SNP_BAITES_FUSA.raw" , skip=2 , na.strings="-")
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
B<-read.table("/NAS/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ/carte_final.txt",header=TRUE)
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
more imputed_hapmap.hmp.txt|sed 's/\<S\>/-/g' |cut -f 1,12- | /NAS/g2pop/HOLTZ_YAN_DATA/programmes/awk_transpose_file.sh | sed 's/ /;/g' | sed 's/rs#/SNP/g'   > fichier_genotypage_QTL.csv


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
scp holtz@CC2-login.cirad.fr://NAS/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/RALIEMENT_BAIT_ET_RNASEQ/IMPUTATION/fichier_genotypage_QTL.csv .


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#











