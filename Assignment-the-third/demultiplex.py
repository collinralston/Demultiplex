#!/usr/bin/env python
import argparse
import bioinfo
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="declaring variables for the file :)")
    parser.add_argument('-1','--filename1',help='your input R1 file',required=True)
    parser.add_argument('-2','--filename2',help='your input R2 file',required=True)
    parser.add_argument('-3','--filename3',help='your input R3 file',required=True)
    parser.add_argument('-4','--filename4',help='your input R4 file',required=True)
    parser.add_argument('-b','--barcodes',help='your input tsv barcodes/indexes file',required=True)
    parser.add_argument('-o','--outpath',help='your output file path',required=True)
    return parser.parse_args()
args = get_args()

indexes = []
with open(args.barcodes,'r') as bc:
    for line in bc:
        line = line.strip()
        if 'sample' in line:
            continue
        else:
            indexes.append(line.split('\t')[4])
barcodes_dict = {}
for i in indexes: # initializing barcodes_dict
    for x in indexes:
        barcodes_dict[f'{i}-{x}'] = 0

#opening all my output files
files = {}
o_path = args.outpath
for i in barcodes_dict:
    b1 = i.split('-')[0]
    b2 = i.split('-')[1]
    if b1 == b2: # Open 48 output fq files in the scratch directory
         files[f'{b1}1']=open(f'{o_path}/{b1}_R1.fq', 'w')
         files[f'{b1}2']=open(f'{o_path}/{b1}_R2.fq', 'w')
files['unk1']=open(f'{o_path}/unknown_R1.fq', 'w')
files['unk2']=open(f'{o_path}/unknown_R2.fq', 'w')
files['hop1']=open(f'{o_path}/hopped_R1.fq', 'w')
files['hop2']=open(f'{o_path}/hopped_R2.fq', 'w')

#opening all my input files
#needs to be gzip.open for final run
in_files = {}
in_files['r1']=open(args.filename1,'rt')
in_files['r2']=open(args.filename2,'rt')
in_files['r3']=open(args.filename3,'rt')
in_files['r4']=open(args.filename4,'rt') # opening all 4 input fq files

def reverse_complement(nuc_str):
    '''returns the reverse complement DNA strand of an input string'''
    comptable = str.maketrans('ACTGN','TGACN')
    comp =  nuc_str.translate(comptable)
    rvs_comp = comp [::-1]
    return rvs_comp
def convert_phred(letter):
    '''Converts a single character into a phred score'''
    return ord(letter) - 33
def print_rec(head,seq,qscore,file):
    '''writes out a single fq record to a specified file'''
    return f'{head}\n{seq}\n+\n{qscore}'
def write_fq_record(rec_list: list,bc1: str,bc2: str,files: dict,file_out: str):
    '''takes the input of all the neccesary components of a fq file
      and rewrites it to a new fq file (already opened)
      bc1 and bc2 are the two barcodes'''
    files[f'{file_out}1'].write(f'{rec_list[0][0]}:{bc1}-{bc2}\n{rec_list[0][1]}\n+\n{rec_list[0][2]}\n')
    files[f'{file_out}2'].write(f'{rec_list[1][0]}:{bc1}-{bc2}\n{rec_list[1][1]}\n+\n{rec_list[1][2]}\n')

unk_count = 0
while True:
    rec_list = [['','',''],['','','']] # list that takes all the record info from R1 and R4 files
    #header lines
    rec_list[0][0]=in_files['r1'].readline().strip()
    if rec_list[0][0]=='':
        break
    rec_list[1][0]=in_files['r4'].readline().strip()
    in_files['r2'].readline()
    in_files['r3'].readline()
    #sequence lines
    seq1 = in_files['r1'].readline().strip()
    rec_list[0][1]=seq1 # saving seq 1 to record
    seq2 = in_files['r2'].readline().strip()
    seq3 = reverse_complement(in_files['r3'].readline()).strip()
    seq4 = in_files['r4'].readline().strip()
    rec_list[1][1]=seq4
    # '+' line
    for i in in_files:
        in_files[i].readline()
    # qscore lines
    q1 = in_files['r1'].readline().strip()
    rec_list[0][2] = q1
    q2 = in_files['r2'].readline().strip()
    q3 = reverse_complement(in_files['r3'].readline()).strip()
    q4 = in_files['r4'].readline().strip()
    rec_list[1][2] = q4
    if seq2 not in indexes or seq3 not in indexes:
        unk_count+=1
        write_fq_record(rec_list,seq2,seq3,files,'unk')
    else:
        badq = False # boolean to track if q score is below threshold
        for i in q2:
            if convert_phred(i) <= 20:
                badq = True
        for i in q3:
            if convert_phred(i) <= 20:
                badq = True
        if badq == True:
            unk_count+=1
            write_fq_record(rec_list,seq2,seq3,files,'unk')
        elif seq2 == seq3:
            barcodes_dict[f'{seq2}-{seq3}']+=1
            write_fq_record(rec_list,seq2,seq3,files,seq2)
        else:
            barcodes_dict[f'{seq2}-{seq3}']+=1
            write_fq_record(rec_list,seq2,seq3,files,'hop')

print('unk count:', unk_count)
hop_count = 0
match_count = 0
for i in barcodes_dict:
    b1 = i.split('-')[0]
    b2 = i.split('-')[1]
    if b1==b2:
        match_count+=barcodes_dict[i]
    else:
        hop_count+=barcodes_dict[i]
print('hop count:',hop_count)
print('match count:', match_count)
with open(f'{o_path}/final_report.txt', 'w') as fr:
    fr.write(f'hop count: {hop_count}')
    fr.write(f'match count: {match_count}')
    fr.write(f'unknown count: {unk_count}')
for i in files:
    files[i].close()
for i in in_files:
    in_files[i].close()

