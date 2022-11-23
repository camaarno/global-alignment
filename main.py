#!/usr/bin/env python

from global_alignment import GlobalAlignment



SEQ_A = "CSTP"
SEQ_B = "CSQTP"

ga = GlobalAlignment()
a, b = ga.align(SEQ_A, SEQ_B)

print(a, b)
