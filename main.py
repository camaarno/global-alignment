#!/usr/bin/env python

from global_alignment import GlobalAligner



SEQ_A = "CQP"
SEQ_B = "CSQATP"

ga = GlobalAligner()
ga.set_blosum('./blosumTest.txt')
result = ga.distance_alignment(SEQ_A, SEQ_B)
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