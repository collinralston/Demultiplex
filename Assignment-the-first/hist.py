#!/usr/bin/env python
import matplotlib.pyplot as plt
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="declaring variables for the file :)")
    parser.add_argument('-f','--filename',help='your input file',required=True)
    parser.add_argument('-o','--output',help='your output file',required=True)
    return parser.parse_args()
args = get_args()
my_list = []
line_ctr = 0
with open(args.filename, 'r') as fh:
    for line in fh:
        line=line.strip()
        my_list.append(float(line))
        line_ctr+=1

        

print(my_list)
x = range(line_ctr)
plt.title('Mean quality score vs. position in sequence')
plt.xlabel('Sequence Position')
plt.ylabel('Mean Phred Quality Score')
plt.bar(x,my_list)
plt.savefig(args.output)