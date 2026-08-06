#!/usr/bin/env python

# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.1"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = 'ACTGN'
RNA_bases = 'ACUGN'

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''

    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Takes a string of phred quality score characters and returns the average quality score"""
    sum=0
    for i in phred_score:
        sum+=convert_phred(i)
    return sum/len(phred_score)

def validate_base_seq(seq: str, RNAFlag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    for i in seq:
        if i not in DNA_bases:
            if i not in RNA_bases:
                return False
    return True

def gc_content(DNA: str):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    if validate_base_seq(DNA) == False:
        return 'invalid sequence'
    DNA = DNA.upper()
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

def calc_median(lst):
    '''Given a sorted list, returns the median value of the list'''
    if len(lst)%2 == 0:
        medians = [lst[len(lst)//2],lst[len(lst)//2-1]]
        return (medians[0]+medians[1])/2
    else:
        return lst[len(lst)//2]

def oneline_fasta(file,out):
    '''takes an input file and an output file - refines the input fasta file for each fasta sequence to only have one line.'''
    with open(file,'r') as file:
        with open(out,'w') as out:
            st = ''
            for line in file:
                if line.startswith('>'):
                    if st != '':
                        out.write(f'{st}\n')
                        st = ''
                    out.write(line)
                else:
                    st += line.strip('\n')
                out.write(st)
            return None


if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    #oneline_test = '>NODE1\nAATGTCGTGC\nAATGG\n>NODE2\nAAAAAAAA\nA'
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    assert calc_median([1,2,3,4,5]) == 3,'median incorrect (odd# of values)'
    assert calc_median([1,2,3,4]) == 2.5,'median incorrect (even# of values)'
    print("Your calc_median function is working! Nice job")
    assert gc_content('AAAGGGCCCTTT') == 0.5, 'GC incorrect'
    assert gc_content('GGGAAATTTA') == 0.3, 'GC incorrect'
    print("Your gc_content function is working! Nice job")
    assert validate_base_seq('eorwiufh') == False, 'this is not a DNA or RNA seq'
    assert validate_base_seq('AAATATGG', False) == True, 'this is a DNA seq'
    assert validate_base_seq('aguauaguc',True) == True, 'this is a RNA seq'
    assert validate_base_seq('AAAATGTGTAC',False) == True, 'this is a DNA seq'
    print("Your validate_base_seq function is working! Nice job")
    assert qual_score('IIIIIIIII') == 40, 'qual score not working - test 1'
    assert qual_score('IC2@') == 30.5, 'qual score not working - test 2'
    print("Your qual_score function is working! Nice job")


