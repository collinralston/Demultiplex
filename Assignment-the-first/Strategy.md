Create dictionary of all index combinations-> index-index - keys, values = # of times that index-index combination occurs.
convert_phred(str)->int: - 
'''takes a  single ascii character and converts it into a quality score'''
    return qscore
    input: I
    output: 40
reverse_complement(str)->str:
'''returns the reverse complement DNA strand of an input string'''
    return reverse_complement
    input: ATGTA
    output:TACAT
Open all R files
    While True loop:
        R4 file, find paired end read
        barcode 1 is R2 file, find barcode
        barcode 2 is the reverse complemented R3 file, find barcode
        if either barcode NOT in the dictionary, if not -> unknown file and +1 count
        else
            if quality score of barcodes is not good enough (convert_phred beforehand for each position in barcode) -> unknown file and count +1
            elif barcodes are equal to each other -> add to that index file count in dict +1
            else - only case left would be if both barcodes are in dict but not matched -> add to barcode combo count in dict -> hopped file
report the number of properly matched indexes, index hopped, and unknown indexes with index_dict counts, hopped count, and unknown count
    to do this loop through index_dict

Defining the problem:
    We are demultiplexing a fastq file, where a bunch of different reads were seqeunced on an Illumina machine, and indexed in order to keep track of whose DNA is whose, and what the source of the DNA is. Our aim is to 'de mix' these reads based on that index, or barcode by creating fq files associated with that index.
Informative output:
    Error rate will be informative. This is expressed in the ratio of the correctly matched indexes vs. the index hopped and unknown indexes.



