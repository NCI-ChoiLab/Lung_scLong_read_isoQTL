#!/bin/sh
# Run TWAS at the isoform level by FUSION
# http://gusevlab.org/projects/fusion/
# Modify from FUSION.compute_weights.R from FUSION
# FILL IN THESE PATHS
module load R
module load plink/1.9.0-beta4.4
module load GCTA
module load gemma
GCTA="/usr/local/apps/GCTA/1.94.3/bin/gcta"
PLINK="/usr/local/apps/plink/1.9.0-beta4.4/plink"
GEMMA="/usr/local/apps/gemma/0.98.5/gemma"
# ALTERNATIVELY: ENSURE THAT plink, gcta, gemma CAN BE CALLED FROM PATH AND REMOVE --PATH_* FLAGS BELOW
# PATH TO DIRECTORY CONTAINING LDREF DATA (FROM FUSION WEBSITE or https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2)
LDREF="/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/bychr"
# THIS IS USED TO RESTRICT INPUT SNPS TO REFERENCE IDS ONLY
while read Cell_type; do
# PATH TO GEUVADIS GENE EXPRESSION MATRIX:
PRE_GEXP=$Cell_type.txt
cd /data/Choi_lung/scLongreads/TWAS/FUSION/$Cell_type/compute_weight
rm -R hsq/
rm -R tmp/
rm -R WEIGHTS/
# GENERATE ID FILE
head -n1 $PRE_GEXP | tr '\t' '\n' | tail -n +5 | awk '{ print $1"\t"$1 }' > ${PRE_GEXP}.ID


# PATH TO OUTPUT DIRECTORY (population-specific subdirs will be made)
OUT_DIR="/data/Choi_lung/scLongreads/TWAS/FUSION/"$Cell_type/"compute_weight/WEIGHTS"
mkdir $OUT_DIR
# ROWS IN THE MATRIX TO ANALYZE (FOR BATCHED RUNS)
BATCH_START=1
# BATCH_END=211

# --- BEGIN SCRIPT:

# NR="${BATCH_START}_${BATCH_END}"
mkdir --parents tmp
mkdir --parents hsq
mkdir --parents out
# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:
mkdir $OUT_DIR

# Loop through each gene expression phenotype in the batch
cat $PRE_GEXP |awk -vs=$BATCH_START 'NR > s' > exp_tmp.txt

while read PARAM; do

# Get the gene positions +/- 1MB
CHR=`echo $PARAM | awk '{ print $1 }'`
P0=`echo $PARAM | awk '{ print $2 - 1e6 }'`
if [ $P0 -lt 0 ]; then
  P0=0
fi
P1=`echo $PARAM | awk '{ print $2 + 1e6 }'`
GNAME=`echo $PARAM | awk '{ print $4 }'`

OUT="$Cell_type.$GNAME"

echo $GNAME $CHR $P0 $P1

# Pull out the current gene expression phenotype
echo $PARAM | tr ' ' '\n' | tail -n +5 | paste $PRE_GEXP.ID - > ./tmp/$OUT.pheno

# Get the locus genotypes for all samples and set current gene expression as the phenotype
$PLINK --bfile $LDREF/genotype.$CHR --pheno ./tmp/$OUT.pheno --make-bed --out ./tmp/$OUT --keep ./tmp/$OUT.pheno --chr $CHR --from-bp $P0 --to-bp $P1 --extract $LDREF/genotype.$CHR.bim

# Process all samples together (for reference purposes only since this is mult-ethnic data)

FINAL_OUT="$OUT_DIR/$GNAME"
cd tmp
Rscript /data/Choi_lung/scLongreads/TWAS/FUSION/fusion_twas-master/FUSION.compute_weights.R --bfile $OUT --tmp $OUT --covar /data/Choi_lung/scLongreads/tensorqtl/isoform_level/$Cell_type/covariates/covs_age_batch_TWAS.txt --out $FINAL_OUT --verbose 0 --hsq_p 0.05 --save_hsq --PATH_plink $PLINK --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet
# ALTERNATIVELY ADD COVARIATES HERE USING THE --covar FLAG
# MINIMAL COMMAND IS: `Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.$pop.tmp --out $FINAL_OUT`
cd ..
# Append heritability output to hsq file
cat $FINAL_OUT.hsq >> hsq/$Cell_type.hsq

# Clean-up just in case
rm -f $FINAL_OUT.hsq $OUT.tmp.*

# Remove all intermediate files
rm $OUT.*

# GO TO THE NEXT GENE
done < "exp_tmp.txt"
done < "celltypes1.txt"

