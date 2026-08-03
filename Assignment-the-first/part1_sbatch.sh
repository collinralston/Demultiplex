#!/usr/bin/env bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
R1=1294_S1_L008_R1_001.fastq.gz
R2=1294_S1_L008_R2_001.fastq.gz
R3=1294_S1_L008_R3_001.fastq.gz
R4=1294_S1_L008_R4_001.fastq.gz
/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/$R3 | grep -A1 '^@' | egrep -v '^@' | grep -v '^--' |grep 'N' | wc -l
#/usr/bin/time -v ./part1.py -f /projects/bgmp/shared/2017_sequencing/$R4 -l 101 -o R4_dist.txt
