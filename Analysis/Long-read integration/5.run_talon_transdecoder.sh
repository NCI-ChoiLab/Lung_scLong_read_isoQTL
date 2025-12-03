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

talon --f config1.csv \
      --db ../myTest32_10s.db \
      --build hg38 \
      --threads 36 \
      --cov 0.95 \
      --identity 0.95 \
      --o output2

module load samtools
module load TransDecoder

pwd

blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep \
	-db ../uniref100/uniprot_sprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
	-num_threads 36 > blastp.outfmt6

set -e
module load hmmer
hmmsearch --cpu $SLURM_CPUS_PER_TASK -E 1e-10 --domtblout pfam.domtblout /data/Choi_lung/scLongreads/TALON_workspace/pfam/Pfam-A.hmm \
	/data/Choi_lung/scLongreads/TALON_workspace/test10s/transcripts.fasta.transdecoder_dir/longest_orfs.pep
