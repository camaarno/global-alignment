#!/usr/bin/env python

from global_alignment import GlobalAlignment



SEQ_A = "CSTP"
SEQ_B = "CSQATP"

ga = GlobalAlignment()
ga.set_blosum('./blosumTest.txt')
ga.distance_alignment(SEQ_A, SEQ_B)


from pprint import pprint

pprint(ga.alignmentMatrix)