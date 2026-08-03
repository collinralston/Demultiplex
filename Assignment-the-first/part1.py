#!/usr/bin/env python
import bioinfo
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="declaring variables for the file :)")
    parser.add_argument('-f','--filename',help='your input file',required=True)
    parser.add_argument('-l','--seqlen',help='your input file sequence length',required=True)
    parser.add_argument('-o','--output',help='your output file',required=True)
    return parser.parse_args()
args = get_args()

def init_list(length: int,lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    i = 0
    if value == []:
        while i<int(length):
            lst.append([])
            i+=1
    else:
        while i < int(length):
            lst.append(value)
            i+=1
    return lst

def populate_list(length,file: str):
    """This function loops through every record in a fastq file and sums all of the quality scores for each base pair. Also keeps a tally for the total
    number of line in the file."""
    with gzip.open(file, 'rt') as fq:
        my_list = init_list(length,[])
        num_lines = 0
        for i,line in enumerate(fq):
            line = line.strip('\n')
            num_lines += 1
            if i%4==3:
                for index,character in enumerate(line):
                    my_list[index]+=bioinfo.convert_phred(character)
            if i%1000000==0:
                print('currently on line', i)
        avg_list = []
        for i in my_list:
            avg_list.append(i/(num_lines/4))
        return avg_list
with open(args.output, 'w') as o:
    lst = populate_list(args.seqlen,args.filename)
    for i in lst:
        o.write(f'{i}\n')


