# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here: [main script](part1.py)     [graphs script](hist.py)        [sbatch script](part1_sbatch.sh)

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | +33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    [R1](plots_and_data/R1_hist.png)
    [R2](plots_and_data/R2_hist.png)
    [R3](plots_and_data/R3_hist.png)
    [R4](plots_and_data/R4_hist.png)
    2. A good quality score cutoff for the index reads and biologcal read pairs could be 20-25. Below the average score in the first base, but not too low to include a bunch of incorrect data. Picking this value now feels like a leap of faith a bit. It may be best to go back and forth, doing multiple runs with different cutoff scores.
    3. R2: 3976613 R3: 3328051
        command: /usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/$R(2 or 3) | grep A1 '^@' | egrep -v '^@' | grep -v '^--' |grep 'N' | wc -l

    
## Part 2
1. Define the problem
    We are demultiplexing a fastq file, where a bunch of different reads were seqeunced on an Illumina machine, and indexed in order to keep track of whose DNA is whose, and what the source of the DNA is. Our aim is to 'de mix' these reads based on that index, or barcode by creating fq files associated with that index.
2. Describe output
    Output will be the 48 r1 and r2 index fq files, as well as the r1 and r2 hopped and unknown fq files. Error rate will be informative. This is expressed in the ratio of the correctly matched indexes vs. the index hopped and unknown indexes.
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode - Done in Strategy.md
5. High level functions. For each function, be sure to include: - Done in Strategy.md [Strategy.md](Mstrategy.md)
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

- Additional Links:
    - [bioinfo module](bioinfo.py)
    - [lab notebook](../notebook_demultiplexBi622.md)

