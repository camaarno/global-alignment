#!/usr/bin/env python

"""
Contains the GlobalDistanceAligner class
"""

"""
@Author: Cameron Arnold
@Data: November 27th, 2022
"""


from operator import lt                      # used to compare scores
from global_aligner_base import GlobalAlignerBase


class GlobalDistanceAligner(GlobalAlignerBase):
    """
    Computes similarity and distance alignments on provided sequences, using
    the affine indel gap model fit sequence alignment
    """

    def __init__(self, blosumPath):
        GlobalAlignerBase.__init__(self, blosumPath)

        self.compare_function = lt
        self.score_function = self.blosum.get_distance_score
        self.gapInitCost = -self.blosum.gapInitCost
        self.gapExtendCost = -self.blosum.gapExtendCost


