- 7/24/26:
    - finished up the unit tests, pseudocode, and answers.md. Will check again in the morning before final submit
- 7/25/26:
    - final upload to github on assignment the first part 1
- 7/28/26:
    - Looked at feedback from my peers on github. I need to do better on the formatting of my pseduocode, and go into more detail generally.
- 7/29/26:
    - Going to work on assignemnt the third in the morning, assignment the first in the afternoon
    - notes on opening files for assignment the third:
        '''
        going to need to open the files outside of the loop, takes a lot more computational power to open a file than to write to it
        '''
    - started demultiplex.py and got all of the output fq files opened and working outside of the main loop. Going to start coding to see if it works on my test files first, then go at it with the real fq files.
    - Going to work on assignment the first part 1: 7:03PM
        Determine which files contain the indexes, and which contain the paired end reads containing the biological data of interest. Create a table and label each file with either read1, read2, index1, or index2.
        | file name | contains |
        | --- | ---|
        |1294_S1_L008_R1_001.fastq.gz|read 1|
        |1294_S1_L008_R2_001.fastq.gz|index 1|
        |1294_S1_L008_R3_001.fastq.gz|index 2|
        |1294_S1_L008_R4_001.fastq.gz|read 2|
        
        command(s) used: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R(1/2/3/4)_001.fastq.gz | head
        note:(1/2/3/4) indicates that I ran the same command for each fq file
        I also determined the which file was which from prior knowledge

        Determine the length of the reads in each file.

        command: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R(1/2/3/4)_001.fastq.gz | head -2 | tail -1 | wc
        result: 102 chars in this sequence line. total sequence is 101bp long (102-newline character)
        note:(1/2/3/4) indicates that I ran the same command for each fq file

        Determine the phred encoding for these data.
    
        command: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -4 | tail -1
        The qscores are phred+33 - The first qscore line contains the char '7', associated with the ascii score 55, which would be impossible with phred+66.
    - realized I didn't have a pixi init when I started working, so I did a pixi init ... I hope this doesn't break anything
    - started the sbatch runs for part1.py ~8:25 PM - going to check in the morning
        R1 jobid# 45868462
            time:46:07.99
            CPU: 99%
            memory:23864
        R2 jobid# 45868463
            time: 8:02.51
            CPU: 99%
            memory: 20936
        R3 jobid# 45868464
            time: 12:22.89
            CPU: 99%
            memory: 20012
        R4 jobid# 45868465
            time: 1:07:56  
            CPU: 99%
            memory: 27040
- 7/30/26
    - runs were successful, it seems like. Got consistent numbers across the board
    - created hist.py to create all the histograms, the graphs look good
- 8/2/26
    - submitted sbatch to see how many indexes have N's
        - command: /usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/$R(2 or 3) | grep A1 '^@' | egrep -v '^@' | grep -v '^--' |grep 'N' | wc -l
        - job 45981926 R2:
            time: 66.86
            CPU: 99%
            memory: 3248
        - job 45981927 R3:
            time: 69.82
            CPU: 99%
            memory: 3236
    
    
    