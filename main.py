#!/usr/bin/env python

import argparse
from os import path
from global_aligner import GlobalAligner
from sequence_reader import SequenceReader


#=========================================================================

# Define our program arguments

def is_valid_file(parser, filePath):
    if path.exists(filePath):
        return filePath
    parser.error("The file %s does not exist!" % filePath)

parser = argparse.ArgumentParser(description='Global Alignment Program')
parser.add_argument('mode', type=str, choices=['distance', 'similarity'],
                    help="The type of global alignment to be performed.")
parser.add_argument('sequences_file', type=lambda x: is_valid_file(parser, x),
                    help="The path to the file containing the two sequences to be compared.")
parser.add_argument('blosum_file', type=lambda x: is_valid_file(parser, x),
                    help="The path to the file containing the BLOSUM matrix and gap penalties.")

args = parser.parse_args()

#=========================================================================

# Run the program

sr = SequenceReader()
seqA, seqB = sr.read_file(args.sequences_file)

ga = GlobalAligner()
ga.set_blosum(args.blosum_file)

result = None
if args.mode == 'similarity':
    alignA, alignB, score = ga.similarity_alignment(seqA, seqB)
else:
    alignA, alignB, score = ga.distance_alignment(seqA, seqB)

print("Sequence A:")
print(alignA, end="\n\n")
print("Sequence B:")
print(alignB, end="\n\n")
print("Score: ", score)




''' TESTING
for row in ga.alignmentMatrix:
    for col in row:
        print(('N/A' if col == None else int(col[0])), "\t", end="")
    print()

print()
print()

for row in ga.deletionMatrix:
    for col in row:
        print(('N/A' if col == None else ('-' if col[1] == None else col[1])), "\t", end="")
    print()
    
'''