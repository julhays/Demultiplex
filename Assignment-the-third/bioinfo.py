#!/usr/bin/env python

# Author: Jules Hays <jkhay@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Calculates the average quality score of a whole FASTQ phred score string"""
    score_sum = 0
    for val in phred_score:
        score = convert_phred(val)
        score_sum = score_sum + score
    return score_sum / len(phred_score)

def validate_base_seq(seq: str, RNAflag: bool = False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(DNA: str) -> float:
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    n = len(lst)
    if n%2 == 1:
        med_position = int((n - 1)/2)
        median = lst[med_position]
    else:
        med_position = int(n/2)
        median = (lst[med_position] + lst[med_position - 1])/2
    return median

def oneline_fasta(filename: str) -> str:
    '''Given a fasta file, returns a new fasta file with each sequence in only 1 line'''
    with open(f'oneln_{filename}', "wt") as fh:
        with open(filename, "rt") as input:
            line1 = True
            for line in input:
                line = line.strip('\n')
                if line1:
                    fh.write(f'{line}\n')
                    line1 = False
                elif line[0] == '>':
                    fh.write(f'\n{line}\n')
                else:
                    fh.write(line)
    return f'oneline_{filename}'

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    #convert_phred tests
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    assert convert_phred("J") == 41, "wrong phred score for 'J'"
    assert convert_phred("u") == 84, "wrong phred score for 'u'"
    assert convert_phred("L") == 43, "wrong phred score for 'L'"
    assert convert_phred("e") == 68, "wrong phred score for 'e'"
    assert convert_phred("S") == 50, "wrong phred score for 'S'"
    print("Your convert_phred function is working! Nice job")

    #qual_score tests
    assert qual_score("EEE") == 36, "wrong average phred score"
    assert qual_score("#I") == 21, "wrong average phred score"
    assert qual_score("EJ") == 38.5, "wrong average phred score"
    assert qual_score("IICC") == 37, "wrong average phred score"
    assert qual_score("JULES") == 44.4, "wrong average phred score"
    assert qual_score("JuLeS") == 57.2, "wrong average phred score"
    print("You calcluated the correct average phred score")

    #valide_base_seq tests
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("JULES", True) == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("AUGTCGTUA") == False, "Validate base seq fails to recognize mix of DNA/RNA"
    assert validate_base_seq("agcttgcnat", False) == True, "Validate base seq dones not work on lowercase DNA"
    print("Passed DNA and RNA tests")

    #gc_content
    assert gc_content("GCGCGC") == 1, "GC contect calculated wrong"
    assert gc_content("AATTATA") == 0, "GC contect calculated wrong"
    assert gc_content("GCATCGAT") == 0.5, "GC contect calculated wrong"
    assert gc_content("GCGCGCTT") == 0.75, "GC contect calculated wrong"
    assert gc_content("GATT") == 0.25, "GC contect calculated wrong"
    assert gc_content("AAAAAAAAAG") == .1, "GC contect calculated wrong"
    print("correctly calculated GC content")

    #calc_median
    assert calc_median([1,2,100]) == 2, "calc_median function does not work for odd length list"
    assert calc_median([1,2]) == 1.5, "calc_median function does not work for even length list"
    assert calc_median([2,7,8,10]) == 7.5, "calc_median function does not work"
    assert calc_median([1,6,200]) == 6, "calc_median function does not work"
    assert calc_median([0,0,0]) == 0, "calc_median function does not work"
    assert calc_median([0.5,2.4,7.8]) == 2.4, "calc_median function does not work for floats"
    print("Median successfully calculated")

    #no assers for oneline_fasta