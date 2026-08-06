#!/usr/bin/env bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
R1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
R2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
R3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
R4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
OUT=/scratch/bgmp/collr/demux/out_files
bcds=/projects/bgmp/shared/2017_sequencing/indexes.txt
/usr/bin/time -v pixi run ./demultiplex.py -1 $R1 -2 $R2 -3 $R3 -4 $R4 -b $bcds -o $OUT -c 25