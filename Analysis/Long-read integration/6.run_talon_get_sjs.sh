#!/bin/sh
#SBATCH --partition=multinode
#SBATCH --constraint=[x2680|x2695]
#SBATCH --nodes=2
#SBATCH --ntasks=36
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mem=108g

##sbatch run_talon.sh

module load samtools
module load mamba_install
source myconda
mamba activate TALON

pwd

talon_get_sjs --gtf /data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoGenes.gtf \
      --ref /data/Choi_lung/lbl/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
      --outprefix NCI_lung_final_sjs
