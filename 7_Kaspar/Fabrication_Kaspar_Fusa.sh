
	#--------------------------------------------------------------------------------------------------
	#   
	#		CREATING KASPAR MARKERS FOR THE FUSA PROJECT
	#
	#					Script réalisé par Yan Holtz (yan1166@hotmail.com / holtz@supagro.inra.fr)
	#
	#---------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------  WHAT IS THIS SCRIPT FOR?

	# This script gives details concerning the 1BS and 5AL QTLs. The output allows the design of corresponding KASPAR Markers.
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#








#------------------------
# 1 --- WHAT AREA DO WE WANT TO TARGET?
# Let's remember the IC we want to target
	#- QTL 1BS: 0-->14.1cM
	#- QTL 5AL: 235.6-->305cM
	
# We use the map with DS bait + RNA-Seq
cd /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/FUSARIOSE/KASPAR
more /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ/map_avec_posi_physique.txt | grep "^1B" | awk '{ if($3 >= 0 && $3<=14.1){print $0,"QTL_1B"}}' > bilan_kaspar_Fusa.txt
more /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/CARTOGRAPHIE_GENETIQUE/2_CARTE_DS/CAPTURE_REF_EPO_ET_RNASEQ/map_avec_posi_physique.txt | grep "^5A" | awk '{ if($3 >= 235.6 && $3<=305){print $0,"QTL_5A"}}' >> bilan_kaspar_Fusa.txt
# I delete 2 SNP that have physical attribution on the 3B chromosome + 1 that makes a bug
more bilan_kaspar_Fusa.txt | grep -v "TRAES3B*" | grep -v "Traes_5AL_1486E0D6C@8" | grep -v "Traes_5AL_F82860F37@31" > tmp
mv tmp bilan_kaspar_Fusa.txt 
# --> 265 markers selected for potential Kaspar creation



#------------------------
# 2 --- ADD CONTIG NAME AND SNP POSITION
more bilan_kaspar_Fusa.txt | awk '{ print $1,$2,$2,$3,$4,$5,$6 }' | sed 's/@/ /' > tmp
mv tmp bilan_kaspar_Fusa.txt
# --> still 265 markers selected for potential Kaspar creation




#------------------------
# 3 -- FIND OUT PARENT'S ALLELES 
python /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/Find_parents_allele_TRAM.py -SNP /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/RNA_SEQ_OCTOBRE_2015/MAPPING_ON_BW/SNP/SNP_super_clean_FUSA_BLE_TENDRE.genot  -out ParAllDS_rnaseq
python /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/Find_parents_allele_TRAM.py -SNP /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/CAPTURE_ALL_INDIV_07_2014/MAPPING_ON_EPO/SNP/SNP_attendus_AND_bonus_clean.genot -out ParAllDS
cat ParAllD* | sort | uniq > ParAllDic
# Can we find difference between genotypes from RNA-Seq and from Capture --> No
more ParAllDic | cut -f1 -d" " | uniq -d
# Merge these allelic parental forms
while read line ; do SNP=$( echo $line | cut -f4 -d" ") ;  a=$(cat ParAllDic | grep -w $SNP | cut -f2,3 -d" ") ;  echo $line" "$a ; done < bilan_kaspar_Fusa.txt > tmp
mv tmp bilan_kaspar_Fusa.txt
# --> still 265 markers selected for potential Kaspar creation
# --> They are located on 138 distinct contigs



#------------------------
# 4 --- GET THE FASTA OF INTERESTING CONTIGS
# list of interesting contigs:
more bilan_kaspar_Fusa.txt | cut -f4 -d" " | sed 's/@.*//' > tmp_liste
# get the fasta sequences:
python /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/take_some_contigs_from_alr_or_fasta.py -fic /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/EPO/0_DATA/RESSOURCES_EPO/EPO_106_After_Homeo_Splitter.fasta -liste tmp_liste -out tmp1
cat /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/TRANSCRIT/RELEASE_22/assemblage_ble_tendre_without_K_D_mart_export.fasta | sed 's/|.*//' > tmp_fasta 
python /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/take_some_contigs_from_alr_or_fasta.py -fic tmp_fasta -liste tmp_liste -out tmp2
cat tmp1 tmp2 > sequence_des_contigs_pour_kaspar.fasta
rm tmp*
# --> 138 Contigs available


#------------------------
# 6 --- IN THIS FASTA, WE WRITE THE SNP BETWEEN BRACKETS
# Get the fasta with [SNP]
more bilan_kaspar_Fusa.txt | cut -f4,9,10 -d" " | sed 's/ AA/ A /g' | sed 's/ CC/ C /g' | sed 's/ GG/ G /g' | sed 's/ TT/ T /g' | tr " " "\t"  > tmp_input_SNP
python /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/make_SNP_file_format_fasta.py -fasta sequence_des_contigs_pour_kaspar.fasta -SNP tmp_input_SNP -format 2 -out tmp_SNP.fasta
more tmp_SNP.fasta_short | grep  -v "^$"  > tmp
mv tmp tmp_SNP.fasta_short
# J'obtiens les fasta complet des contigs avec les SNPs dedans + des courts fasta autours du SNP
# Ajout des fasta au bilan + Ajout en tete!
while read line ; do a=$(echo $line| cut -f2 -d" "  ) ; b=$(cat tmp_SNP.fasta | grep -A1 $a | tail -1) ; echo $line $b  ; done < bilan_kaspar_Fusa.txt > tmp
# Idem pour les fasta short
while read line ; do a=$(echo $line| cut -f4 -d" "  ) ; b=$(cat tmp_SNP.fasta_short | grep -A1 $a | tail -1) ; echo $line $b  ; done < tmp > tmp1
mv tmp1 bilan_kaspar_Fusa.txt
rm ParAllD*



#------------------------
# 7 --- BLAST SEQUENCES TO DETECT HOMEOLOGOUS SEQUENCES
#On va blaster les seq entournt les SNP contenant les marqueurs QTL sur notre référence fasta, pour voir si il y a des blasts très similaires!
#Blast sur le BW
more tmp_SNP.fasta_short | sed 's/\[//' | sed 's/\/.\]//' > tmp_fasta_short
makeblastdb -in tmp_fasta_short -dbtype nucl -out database
blastn -db database  -query  /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/TRANSCRIT/RELEASE_28/transcriptome_BW_without_K_D.fasta -outfmt '6 qseqid sseqid qstart qend qlen sstart send slen length pident' > resultat_blastn
#Filtre a au moins 90% de similarité
more resultat_blastn | awk '{ if($10>90){print $0}}' > clean_blast_result
rm database* resul*

# Ajout des résultat de blast au bilan
while read line ; do a=$(echo $line| cut -f4 -d" " ) ; b=$(cat clean_blast_result | grep -w $a | grep "^Traes_[1-7]A[L,S]_" | cut -f1,10 | sort | tail -1)  ; [  -z "$b" ] && b=$(echo "- -")  ;   c=$(cat clean_blast_result | grep $a | grep "^Traes_[1-7]B[L,S]_"| cut -f1,10 | sort | tail -1)  ;   [  -z "$c" ] && c=$(echo "- -") ; echo $line "@@"$b "@@"$c  ; done < bilan_kaspar_Fusa.txt > tmp



# 8 --- RECUPERATION DES SEQ DE BLAST
more clean_blast_result | cut -f1 > liste
python /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/take_some_contigs_from_alr_or_fasta.py  -fic /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/BREAD_WHEAT_IWGSC_and_HORDEUM/TRANSCRIT/RELEASE_28/transcriptome_BW_without_K_D.fasta -liste liste -out contigs_with_blast.fasta


# 9 --- ADD HEADER + CUSTOM FORMAT
echo "Chromosome Contig Position Marker Position_cM Chromo_pos_physique Position_BP Nature Allele_parent1 Allele_Parent2 contig_sequence_long contig_sequence_short blast_genomeA similarity_genomeA blast_genomeB similarity_genomeB" > bilan_kaspar_Fusa.txt
cat tmp | sed 's/@@ /- /g' | sed 's/@@//g' | sort -k5n -t" " >> bilan_kaspar_Fusa.txt
more bilan_kaspar_Fusa.txt    | sed 's/ /;/g' > bilan_kaspar_Fusa.csv
rm tmp* clean_blast_result liste bilan_kaspar_Fusa.txt
# Transfert en local
cd /Users/holtz/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/SCRIPT/6_Kaspar
scp holtz@CC2-login.cirad.fr:/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/FUSARIOSE/KASPAR/bila* .
scp holtz@CC2-login.cirad.fr:/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/FUSARIOSE/KASPAR/contigs_with*fasta .

# Caractérisation vite fait
R
data=read.table( "bilan_kaspar_Fusa.csv" , sep=";" , header=T)
#Nbr de Kaspar possible:
nrow(data)
#Nbr pour le 1B et pour le 1A:
table(data$Chromosome)
#Nbr de position unique pour le 1B et pour le 1A
length(unique(data[data$Chromosome=="1B",5]))
length(unique(data[data$Chromosome=="5A",5]))



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


