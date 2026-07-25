# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |  |  |
| 1294_S1_L008_R2_001.fastq.gz |  |  |  |
| 1294_S1_L008_R3_001.fastq.gz |  |  |  |
| 1294_S1_L008_R4_001.fastq.gz |  |  |  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem
    We are demultiplexing a fastq file, where a bunch of different reads were seqeunced on an Illumina machine, and indexed in order to keep track of whose DNA is whose, and what the source of the DNA is. Our aim is to 'de mix' these reads based on that index, or barcode by creating fq files associated with that index.
2. Describe output
    Output will be the 48 r1 and r2 index fq files, as well as the r1 and r2 hopped and unknown fq files. Error rate will be informative. This is expressed in the ratio of the correctly matched indexes vs. the index hopped and unknown indexes.
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode - Done in Strategy.md
5. High level functions. For each function, be sure to include: - Done in Strategy.md
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
