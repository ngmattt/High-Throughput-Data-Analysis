#!/bin/bash
#SBATCH -J oases_multiK
#SBATCH -A r00270
#SBATCH -p gpu
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 18:00:00
#SBATCH -o oases_multiK.%j.out
#SBATCH -e oases_multiK.%j.err


export PS1=${PS1:-"$"}
set +u
module load conda
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate assignment1

# Parameters
READ1="/N/u/ngmat/Quartz/assignment1/trimming/SRR21904868_1_val_1.fq"
READ2="/N/u/ngmat/Quartz/assignment1/trimming/SRR21904868_2_val_2.fq"
KLIST="41 51 61 71 81 91 101"
INS=300
THREADS="${SLURM_CPUS_PER_TASK:-12}"

########################
: > oases_results.txt  # initialize summary file once

for K in $KLIST; do
  OUT="oases_${K}"
  echo "Running Oases transcript assembler with k-mer: $K"
  mkdir -p "$OUT"

  velveth "$OUT" "$K" -shortPaired -fastq -separate "$READ1" "$READ2"
  velvetg "$OUT" -read_trkg yes
  oases   "$OUT" -ins_length "$INS" -min_trans_lgth 200 -cov_cutoff 3

  # Summary Statistics
  F="$OUT/transcripts.fa"
  echo "K-mer $K results:" >> oases_results.txt
  echo "Number of transcripts: $(grep -c '^>' "$F")" >> oases_results.txt
  echo "Total length (bp): $(grep -v '^>' "$F" | tr -d '\n' | wc -c)" >> oases_results.txt
  n50=$(awk '/^>/{if(L){print L};L=0;next}{L+=length}END{if(L)print L}' "$F" \
        | sort -nr \
        | awk '{a[NR]=$1;s+=$1} END{h=s/2;c=0; for(i=1;i<=NR;i++){c+=a[i]; if(c>=h){print a[i]; exit}}}')
  echo "N50: ${n50:-0}" >> oases_results.txt
  echo "---" >> oases_results.txt

  echo "[k=$K] transcripts: $F"
done
