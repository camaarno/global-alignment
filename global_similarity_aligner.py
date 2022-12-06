#!/usr/bin/env python

"""
Contains the GlobalSimilarityAligner class
"""

"""
@Author: Cameron Arnold
@Data: November 27th, 2022
"""


from operator import gt                      # used to compare scores
from global_aligner_base import GlobalAlignerBase


class GlobalSimilarityAligner(GlobalAlignerBase):
    """
    Computes similarity alignments on provided sequences, using
    the affine indel gap model fit sequence alignment
    """

    def __init__(self, blosumPath):
        GlobalAlignerBase.__init__(self, blosumPath)

        self.compare_function = gt
        self.score_function = self.blosum.get_similarity_score
        self.gapInitCost = self.blosum.gapInitCost
        self.gapExtendCost = self.blosum.gapExtendCost

