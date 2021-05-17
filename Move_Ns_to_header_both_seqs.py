import numpy as np

#this function to look for the distance between two strings was taken from:
#https://www.datacamp.com/community/tutorials/fuzzy-string-python
def levenshtein_ratio_and_distance(s, t, ratio_calc = False):
    """ levenshtein_ratio_and_distance:
        Calculates levenshtein distance between two strings.
        If ratio_calc = True, the function computes the
        levenshtein distance ratio of similarity between two strings
        For all i and j, distance[i,j] will contain the Levenshtein
        distance between the first i characters of s and the
        first j characters of t
    """
    # Initialize matrix of zeros
    rows = len(s)+1
    cols = len(t)+1
    distance = np.zeros((rows,cols),dtype = int)

    # Populate matrix of zeros with the indeces of each character of both strings
    for i in range(1, rows):
        for k in range(1,cols):
            distance[i][0] = i
            distance[0][k] = k

    # Iterate over the matrix to compute the cost of deletions,insertions and/or substitutions
    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0 # If the characters are the same in the two strings in a given position [i,j] then the cost is 0
            else:
                # In order to align the results with those of the Python Levenshtein package, if we choose to calculate the ratio
                # the cost of a substitution is 2. If we calculate just distance, then the cost of a substitution is 1.
                if ratio_calc == True:
                    cost = 2
                else:
                    cost = 1
            distance[row][col] = min(distance[row-1][col] + 1,      # Cost of deletions
                                 distance[row][col-1] + 1,          # Cost of insertions
                                 distance[row-1][col-1] + cost)     # Cost of substitutions
    if ratio_calc == True:
        # Computation of the Levenshtein Distance Ratio
        Ratio = ((len(s)+len(t)) - distance[row][col]) / (len(s)+len(t))
        return Ratio
    else:
        # print(distance) # Uncomment if you want to see the matrix showing how the algorithm computes the cost of deletions,
        # insertions and/or substitutions
        # This is the minimum number of edits needed to convert string a to string b
        #return "The strings are {} edits away".format(distance[row][col])
        return distance[row][col]

from sys import argv
import re #regular expressions package
import gzip
from itertools import islice
#right now this takes about 30 minutes to run for a 3GB fastq.gz file.
#The line counting alone takes 8 minutes. Try to parallelize this?
#this works for Python version 2.7 there are some new features in Python 3 that
    #will break parts of this script

script, InFileName, InFileName2 = argv

OutFileName = InFileName + 'edited'
OutFileName2 = InFileName2 + 'edited'


def file_len(fname):
    with gzip.open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
count2 = file_len(InFileName)
#print(count2)
seq = count2/4
print("There are %s sequences.") % seq


iterations = 0
eleven_count = 0
Ncount = []
seqNames = []
Ns = []


def next_n_lines(file_opened, N):
    return [x.strip() for x in islice(file_opened, N)]

#OutFile = open(OutFileName, 'w')
with open(OutFileName, 'w') as OutFile, gzip.open(InFileName, 'r') as infile:
    while iterations < seq:
        lines = next_n_lines(infile, 4)
        #print(lines[0].strip('\n'))
        if "N" in lines[1][0:5]:
            Ncount.append(iterations)
            seqNames.append(lines[0])

        # elif lines[1][6:10]=='CACG':
        #     Line = lines[0].strip('\n')
        #     Line = re.sub(r'\:CATGCCTA\+.+', '', Line)
        #     NewLine = (Line + ':' + lines[1][0:6] + '\n' )
        #     OutFile.write(NewLine)
        #     OutFile.write(lines[1][6:]+'\n')
        #     OutFile.write(lines[2]+'\n')
        #     OutFile.write(lines[3][6:]+'\n')
        #     eleven_count += 1
        else:
            Line = lines[0].strip('\n')
            Line = re.sub(r' 1\:N.+', '', Line)
            NewLine = (Line + ' 1:N:0:' + lines[1][0:5] + '\n' )
            Ns.append(lines[1][0:5])
            OutFile.write(NewLine)
            OutFile.write(lines[1][5:]+'\n')
            OutFile.write(lines[2]+'\n')
            OutFile.write(lines[3][5:]+'\n')

        iterations += 1
print("There were %s rescues for barcode 11.") % eleven_count
print("There were %s Ns in unique identifiers") % len(Ncount)
percent=(len(Ncount)*100/seq)
print("%s percent have Ns in the unique identifiers") % percent
#print(Ncount)
#print(seqNames)

count3 = file_len(InFileName2)
#print(count3)
print("check that the file lengths match:")
print(count2 == count3)

iterations2 = 0
remove = 0
test = 0
with open(OutFileName2, 'w') as OutFile, gzip.open(InFileName2, 'r') as infile2:
    while iterations2 < seq:
        if iterations2 in Ncount:
            lines = next_n_lines(infile2, 4)
            if levenshtein_ratio_and_distance(seqNames[remove], lines[0])>=2:
                print('files do not match!')
                break
            else:
#                print('removing seq')
                remove += 1
        else:
            lines = next_n_lines(infile2, 4)
            Line = lines[0].strip('\n')
            Line = re.sub(r' 2\:N.+', '', Line)
            NewLine = (Line + ' 2:N:0:' + Ns[test] + '\n' )
            OutFile.write(NewLine)
            OutFile.write(lines[1]+'\n')
            OutFile.write(lines[2]+'\n')
            OutFile.write(lines[3]+'\n')
            test +=1

        iterations2 += 1
