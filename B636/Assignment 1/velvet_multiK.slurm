#!/bin/bash
#SBATCH -J velvet_multiK
#SBATCH -A r00270
#SBATCH -p gpu
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 18:00:00
#SBATCH -o velvet_multiK.%j.out
#SBATCH -e velvet_multiK.%j.err

export PS1=${PS1:-"$"}
set +u
module load conda
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate assignment1
########################
# Parameters
READ1="/N/u/ngmat/Quartz/assignment1/trimming/SRR21904868_1_val_1.fq"
READ2="/N/u/ngmat/Quartz/assignment1/trimming/SRR21904868_2_val_2.fq"
KLIST="41 51 61 71 81 91 101"
INS=300
INSSD=50
THREADS="${SLURM_CPUS_PER_TASK:-12}"

########################
: > velvet_results.txt

for K in $KLIST; do
  OUT="velvet_${K}"
  echo "Running velvet assembler with k-mers: $KLIST"
  mkdir -p "$OUT"

  # Build hash
  velveth "$OUT" "$K" -shortPaired -fastq -separate "$READ1" "$READ2"

  # Graph + assembly
  velvetg  "$OUT" \
    -exp_cov auto \
    -cov_cutoff auto \
    -ins_length "$INS" \
    -ins_length_sd "$INSSD" \
    -scaffolding yes \
    -min_contig_lgth 0 \
    -read_trkg yes

# Summary Statistics
  echo "K-mer $K results:" >> velvet_results.txt
  echo "Number of contigs: $(grep -c '^>' "$OUT/contigs.fa")" >> velvet_results.txt
  echo "Total length (bp): $(grep -v '^>' "$OUT/contigs.fa" | tr -d '\n' | wc -c)" >> velvet_results.txt
  n50=$(awk '/^>/{if(L){print L};L=0;next}{L+=length}END{if(L)print L}' "$OUT/contigs.fa" \
        | sort -nr \
        | awk '{a[NR]=$1;s+=$1} END{h=s/2;c=0; for(i=1;i<=NR;i++){c+=a[i]; if(c>=h){print a[i]; exit}}}')
  echo "N50: ${n50:-0}" >> velvet_results.txt
  echo "---" >> velvet_results.txt
done
