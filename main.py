#!/usr/bin/env python

import argparse
from os import path
from sequence_reader import SequenceReader
from global_aligner import GlobalAligner


MATCH_SYMBOL = '*'
MISMATCH_SYMBOL = '|'
GAP_SYMBOL = ' '

#=========================================================================

# Define our program arguments

def is_valid_file(parser, filePath):
    if path.exists(filePath):
        return filePath
    parser.error("The file %s does not exist!" % filePath)

parser = argparse.ArgumentParser(description='Global Alignment Program')
parser.add_argument('mode', type=str, choices=['distance', 'similarity'],
                    help="The type of global alignment to be performed.")
parser.add_argument('blosum_file', type=lambda x: is_valid_file(parser, x),
                    help="The path to the file containing the BLOSUM matrix and gap penalties.")
parser.add_argument('sequences_file', type=lambda x: is_valid_file(parser, x),
                    help="The path to the file containing the two sequences to be compared.")

args = parser.parse_args()

#=========================================================================

# Run the program

sr = SequenceReader()
sr.set_sequences(args.sequences_file)

ga = GlobalAligner()
ga.set_blosum(args.blosum_file)

alignment_function = None
if args.mode == 'similarity':
    alignment_function = ga.similarity_alignment
else:
    alignment_function = ga.distance_alignment

for i in range(sr.count - 1):
    for j in range(i + 1, sr.count):
        seqA = sr.sequences[i]
        seqB = sr.sequences[j]
        alignA, alignB, score = alignment_function(seqA, seqB)

        # Print names
        print(f"Alignment of {sr.proteinNames[i]} ({sr.speciesNames[i]}) and {sr.proteinNames[j]} ({sr.speciesNames[j]})")

        # Print the first alignment
        print(alignA)

        # Print a row of symbols representing character matches
        symbol = ''
        for k in range(len(alignA)):
            if alignA[k] == alignB[k]:
                symbol = MATCH_SYMBOL
            elif alignA[k] == ga.GAP_SYMBOL or alignB[k] == ga.GAP_SYMBOL:
                symbol = GAP_SYMBOL
            else:
                symbol = MISMATCH_SYMBOL
            print(symbol, end="")
        print()

        # Print the second alignment
        print(alignB)

        # Print the score
        print("Score: ", score)
        print()


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