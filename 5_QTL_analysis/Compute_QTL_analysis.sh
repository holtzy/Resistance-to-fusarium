	#--------------------------------------------------------------------------------------------------
	#   
	#		QTL detection for the fusariose resistance
	#
	#					Script réalisé par Yan Holtz (yan1166@hotmail.com / holtz@supagro.inra.fr)
	#
	#---------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------  SCRIPT GOAL

# This script is made to compute QTL using QTLREL

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# QTL FUSA

# COPY GENO, PHENO and MAP on CC2:
cd /Users/yan/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA
scp fichier_genotypage_QTL.csv map_avec_posi_physique.txt	PHENOTYPE/phenotypage_all_fusa_blup.csv holtz@CC2-login.cirad.fr:/homedir/holtz/work/FUSA/QTL


# Compute QTL with QTLREL in CC2. Envoie en 10 fois?
cd  /homedir/holtz/work/FUSA/QTL
# a partir de quel colonne je travaille
a=2
# Combien de caractère par qsub
b=5
for i in $(seq 1 36); do
	mkdir TMP_PART$i
	cat phenotypage_all_fusa_blup.csv | cut -f1,$a-$b -d";" > TMP_PART$i/fic_pheno.csv
	cd TMP_PART$i
	qsub -q normal.q -b yes -cwd -N tmp_qtl "Rscript /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/QTL_detection_with_QTLRel.R  ../fichier_genotypage_QTL.csv fic_pheno.csv  ../map_avec_posi_physique.txt"
	cd ..
	((a=$a+14))
	((b=$b+14))
done

# Check que tout a bien marché:
for i in TMP* ; do echo $i ; ls $i/bilan* ; done

# Regroupe les fichiers
more TMP_PART1/bilan_simple_marker | head -1 > bilan_simple_marker
for i in TMP_*/bila* ; do cat $i | sed '1d' >> bilan_simple_marker ; done
gzip bilan_simple_marker

# Transfert results in the dropbox
cd ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA
scp holtz@CC2-login.cirad.fr:~/work/FUSA/QTL/bilan_simple_marker.gz .

# --> and use the RMD file to analyse these QTLs!

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# QTL METABOLOMIQUE

# COPY GENO, PHENO and MAP on CC2:
cd /Users/yan/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA
scp fichier_genotypage_QTL.csv map_avec_posi_physique.txt	PHENOTYPE/phenotypage_all_metabolomique.csv  CC2-login.cirad.fr:/gs7k1/projects/g2pop/HOLTZ_YAN_DATA/FUSARIOSE/QTL_METABOLOMIQUE

# Compute QTL with QTLREL in CC2. Envoie en 10 fois?
cd  /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/FUSARIOSE/QTL_METABOLOMIQUE
a=2
b=15
for i in $(seq 1 20); do
	mkdir TMP_PART$i
	cat phenotypage_all_metabolomique.csv | cut -f1,$a-$b -d";" > TMP_PART$i/fic_pheno.csv
	cd TMP_PART$i
	qsub -q normal.q -b yes -cwd -N tmp_qtl "Rscript /gs7k1/projects/g2pop/HOLTZ_YAN_DATA/programmes/QTL_detection_with_QTLRel.R  ../fichier_genotypage_QTL.csv fic_pheno.csv  ../map_avec_posi_physique.txt"
	cd ..
	((a=$a+14))
	((b=$b+14))
done

# Regroupe les fichiers
more TMP_PART1/bilan_simple_marker | head -1 > bilan_simple_marker
for i in TMP_*/bila* ; do cat $i | sed '1d' >> bilan_simple_marker ; done

# Transfert results in the dropbox
cd ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA
scp holtz@CC2-login.cirad.fr://gs7k1/projects/g2pop/HOLTZ_YAN_DATA/FUSARIOSE/QTL_METABOLOMIQUE/bilan_simple_marker bilan_simple_marker_metabolomique


# --> and use the RMD file to analyse these QTLs!

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#







#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Use the appli quickly
cd ~/Dropbox/APPLI_SHINY_QTL/DATA_FUSA
cp ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA/bilan_simple_marker .
cp ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA/map_avec_posi_physique.txt carte
cp ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA/fichier_genotypage_QTL.csv genotypage.csv
cp ~/Dropbox/Publi_Fusariose/ANALYSIS_REPRO/DATA/PHENOTYPE/phenotypage_all_fusa.csv phenotypage.csv


cd /Users/yan/Dropbox/APPLI_SHINY_QTL
R
library(shiny)
runApp("SHINY_APP_FOR_QTL_ANALYSIS")
















