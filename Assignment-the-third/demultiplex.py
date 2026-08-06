#!/usr/bin/env python
import argparse
import bioinfo
import gzip
import matplotlib.pyplot as plt
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="declaring variables for the file :)")
    parser.add_argument('-1','--filename1',help='your input R1 file',required=True)
    parser.add_argument('-2','--filename2',help='your input R2 file',required=True)
    parser.add_argument('-3','--filename3',help='your input R3 file',required=True)
    parser.add_argument('-4','--filename4',help='your input R4 file',required=True)
    parser.add_argument('-b','--barcodes',help='your input tsv barcodes/indexes file',required=True)
    parser.add_argument('-o','--outpath',help='your output file path',required=True)
    parser.add_argument('-c','--cutoff',help='your phred score cutoff',required=True)
    return parser.parse_args()
args = get_args()
cutoff = int(args.cutoff)

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
in_files['r1']=gzip.open(args.filename1,'rt')
in_files['r2']=gzip.open(args.filename2,'rt')
in_files['r3']=gzip.open(args.filename3,'rt')
in_files['r4']=gzip.open(args.filename4,'rt') # opening all 4 input fq files

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
            if convert_phred(i) <= cutoff:
                badq = True
        for i in q3:
            if convert_phred(i) <= cutoff:
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

array = np.zeros((24,24),dtype=float)
row=0
index=[]
col=0
ct = 0
sum=hop_count+match_count+unk_count
match_dict = {}
with open('summary.md','w') as sm:
    sm.write(f'unknown index count: {unk_count}')
    sm.write(f'hopped index count: {hop_count}')
    sm.write(f'matched index count: {match_count}')
    for key in barcodes_dict:
        b1 = key.split('-')[0]
        b2 = key.split('-')[1]
        if b1==b2:
            sm.write(f'{b1} percentage: {(barcodes_dict[key]/sum)*100}\n')
            match_dict[b1]= barcodes_dict[key]
            print(f'{b1} percentage: {(barcodes_dict[key]/sum)*100}')
        if col == 24:
            row+=1
            col=0
        if col == 0:
            index.append(key.split('-')[0]) # y coordinate
        array[row,col] = np.log((barcodes_dict[key]+1))
        #print('row',row)
        #print('col',col)
        col+=1
    sm.write(f'total match percentage: {(match_count / sum)*100}\n')
    sm.write(f'hopping percentage: {(hop_count / sum)*100}\n')
    sm.write(f'unknown percentage: {(unk_count / sum)*100}')
#https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
fig,ax = plt.subplots(figsize=(10.,10.))
im=ax.imshow(array,cmap='cividis') #https://stackoverflow.com/questions/32236046/add-a-legend-to-my-heatmap-plot
plt.colorbar(im,label='ln(#of index pairs +1)')
ax.set_xticks(range(len(index)),labels=index,
              rotation=45, rotation_mode="xtick")
ax.set_yticks(range(len(index)),labels=index)
ax.set_title('Barcode Matching and Hopping in a Illumina Sequencing Run')
ax.set_xlabel('Index 2')
ax.set_ylabel('Index 1')
plt.savefig('heatmap.png')
sum=hop_count+match_count+unk_count

hop_dict = {}
for key in barcodes_dict:
    b1 = key.split('-')[0]
    b2 = key.split('-')[1]
    if b1!=b2:
        if b2 not in hop_dict and b1 not in hop_dict:
            hop_dict[b1] = barcodes_dict[key]
            hop_dict[b2] = barcodes_dict[key]
        elif b2 in hop_dict and b1 not in hop_dict:
            hop_dict[b1] = barcodes_dict[key]
            hop_dict[b2] += barcodes_dict[key]
        elif b1 in hop_dict and b2 not in hop_dict:
            hop_dict[b2] = barcodes_dict[key]
            hop_dict[b1] += barcodes_dict[key]
        else:
            hop_dict[b1] += barcodes_dict[key]
            hop_dict[b2] += barcodes_dict[key]
diff_dict={}
for i in hop_dict:
    for j in match_dict:
        if j==i:
            diff_dict[i]=np.abs(match_dict[i]/match_count*100 - hop_dict[i]/hop_count*100)
diff_dict = dict(sorted(diff_dict.items(),key=lambda item: item[1],reverse=True))
fig,ax = plt.subplots(figsize = (10.,10.))
ax.set_title('Difference in Hopping and Matching per Index')
ax.set_ylabel('difference in percentage between hopping and matching')
ax.set_xlabel('index')
ax.set_xticks(range(len(diff_dict)),labels=list(diff_dict.keys()),
              rotation=45, rotation_mode="xtick")
plt.bar(list(diff_dict.keys()),list(diff_dict.values()))
plt.savefig('diff_HopAndMatch.png')

#this prints out two plots, side by side, that compare the matched vs hopped barcodes to see
# a side by side comparison, similar to diff.png but without the datasets combined.
#https://stackoverflow.com/questions/613183/how-do-i-sort-a-dictionary-by-value
#https://stackoverflow.com/questions/22799308/how-can-i-save-two-plots-on-a-single-file-in-python
#hop_dict = dict(sorted(hop_dict.items(),key=lambda item: item[1],reverse=True))
#match_dict = dict(sorted(match_dict.items(),key=lambda item: item[1],reverse=True))
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(20.,20.))
ax1.set_ylabel('hopping count')
ax1.set_xlabel('index')
ax1.set_xticks(range(len(hop_dict)),labels=list(hop_dict.keys()),
              rotation=45, rotation_mode="xtick")
ax1.set_title('Index hopping in a multiplexed Illumina run')
ax1.bar(list(hop_dict.keys()),list(hop_dict.values()))
ax2.set_ylabel('matching count')
ax2.set_xlabel('index')
ax2.set_xticks(range(len(match_dict)),labels=list(match_dict.keys()),
              rotation=45, rotation_mode="xtick")
ax2.set_title('Index matching in a multiplexed Illumina run')
ax2.bar(list(match_dict.keys()),list(match_dict.values()))
plt.savefig('hops_and_matches.png')







    
    
        
    
    


