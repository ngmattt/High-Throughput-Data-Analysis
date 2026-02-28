The entire assignment was done in Quartz hpc. 
The following modules were loaded for the assignment:
conda, python, sra-toolkit, trimgalore, cutadapt, velvet, oases 

Conda environment was created for this assignment:
    'conda create --name assignment1'
    'conda activate assignment1'

To access the data, I used the following commands: 
    'fasterq-dump --split-files SRR21904868'  

Assembly tools were installed via bioconda into the environment:
    'conda install bioconda::sra-tools velvet oases'

For preprocessing, fastqc and trimgalore was used to evaluate and clean the fastq files:
    'fastqc -o fastqc SRR21904868_1.fastq SRR21904868_2.fastq'
    'mkdir fastqc'
    'fastqc -o fastqc SRR21904868_1.fastq SRR21904868_2.fastq'
    'trim_galore --phred33 --fastqc --paired SRR21904868_1.fastq SRR21904868_2.fastq -o trimming'

For the genome assembly, two slurm job scripts were created for velvet and oases. 

The following code block is the start of the slurm scripts which is used to request a predetermined amount of resources to run the jobs. It also includes parameters used by the tools to perform the assembly such as the preprocessed fastq files and INS & INSSD which are metrics that indicate read length and standard deviation. 

```bash
#!/bin/bash
#SBATCH -J velvet_multiK
#SBATCH -A r00270
#SBATCH -p gpu
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 18:00:00
#SBATCH -o velvet_multiK.%j.out
#SBATCH -e velvet_multiK.%j.err

####################################################################################
# Loading the environment and preparing the parameters for the genome assembly task
export PS1=${PS1:-"$"}
set +u
module load conda
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate assignment1

READ1="/N/u/ngmat/Quartz/assignment1/trimming/SRR21904868_1_val_1.fq"
READ2="/N/u/ngmat/Quartz/assignment1/trimming/SRR21904868_2_val_2.fq"
KLIST="41 51 61 71 81 91 101"
INS=300
INSSD=50
THREADS="${SLURM_CPUS_PER_TASK:-12}"
```

The following code block shows the commands ran for the velvet tool and output with summary statistics including K-mer length, number of contigs, and total length (bp) into a .txt file called 'velvet_results.txt'.

```bash
: > velvet_results.txt # initialize the results summary file

for K in $KLIST; do # iterate over the list of k-mers
  OUT="velvet_${K}"
  echo "Running velvet assembler with k-mers: $KLIST"
  mkdir -p "$OUT"

  # Velvet tool 1st stage (fastq reads as input - hash table used to construct graph as output)
  velveth "$OUT" "$K" -shortPaired -fastq -separate "$READ1" "$READ2" 

  # Velvet tool 2nd stage (graph & assembly)
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
```
The following code block shows the commands ran for the oases tool and output with summary statistics including K-mer length, number of transcripts, and total length (bp) into a .txt file called 'oases_results.txt'.

```bash
: > oases_results.txt  # initialize results summary file

for K in $KLIST; do
  OUT="oases_${K}"
  echo "Running Oases transcript assembler with k-mer: $K"
  mkdir -p "$OUT"

  velveth "$OUT" "$K" -shortPaired -fastq -separate "$READ1" "$READ2" 
  velvetg "$OUT" -read_trkg yes
  oases   "$OUT" -ins_length "$INS" -min_trans_lgth 200 -cov_cutoff 3 # Run Oases on Velvet’s graph to assemble transcripts

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
```

Results:

The slurm scripts generated a summary .txt file that includes the number of contigs/transcripts, total length (bp) of assembly, and N50 values for each k-mer length. 

velvet_results.txt
```txt
K-mer 41 results:
Number of contigs: 316
Total length (bp): 4756813
N50: 67950
---
K-mer 51 results:
Number of contigs: 261
Total length (bp): 4763624
N50: 92070
---
K-mer 61 results:
Number of contigs: 255
Total length (bp): 4771103
N50: 94426
---
K-mer 71 results:
Number of contigs: 256
Total length (bp): 4779838
N50: 112548
---
K-mer 81 results:
Number of contigs: 239
Total length (bp): 4784685
N50: 118766
---
K-mer 91 results:
Number of contigs: 227
Total length (bp): 4786706
N50: 112588
---
K-mer 101 results:
Number of contigs: 285
Total length (bp): 4804038
N50: 101753
---
```

oases_results.txt
```txt
K-mer 41 results:
Number of transcripts: 502
Total length (bp): 10618537
N50: 54116
---
K-mer 51 results:
Number of transcripts: 388
Total length (bp): 10629778
N50: 84737
---
K-mer 61 results:
Number of transcripts: 404
Total length (bp): 12121704
N50: 143928
---
K-mer 71 results:
Number of transcripts: 400
Total length (bp): 11852008
N50: 310285
---
K-mer 81 results:
Number of transcripts: 442
Total length (bp): 10631981
N50: 143968
---
K-mer 91 results:
Number of transcripts: 367
Total length (bp): 12169471
N50: 553724
---
K-mer 101 results:
Number of transcripts: 595
Total length (bp): 12636781
N50: 257507
---
```
