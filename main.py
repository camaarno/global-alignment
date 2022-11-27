#!/usr/bin/env python

from global_aligner import GlobalAligner
from sequence_reader import SequenceReader



sr = SequenceReader()
seqA, seqB = sr.read_file('./sequences.txt')

ga = GlobalAligner()
ga.set_blosum('./blosumTest.txt')
result = ga.distance_alignment(seqA, seqB)
print(result)

'''
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