#!/usr/bin/env python

from blosum import Blosum

from pprint import pprint


class GlobalAlignment:


    def __init__(self):
        self.alignmentMatrix = []   # optimal alignment matrix     (A)
        self.deletionMatrix = []    # alignment ends with deletion (D)
        self.insertionMatrix = []   # alignment ends with deletion (I)
        self.seqA = ""              # horizontal sequence
        self.seqB = ""              # vertical sequence

        self.matrixHeight = 0   # height of each matrix
        self.matrixWidth = 0    # width of each matrix

        self.blosum = Blosum()
        self.gapInitCost = 0    # gap initiation cost
        self.gapExtendCost = 0  # gap extension cost

        self.func = None  # min / max
        self.scorer = None  # blosum.get_similarity_score / blosum.get_distance_score


    def set_blosum(self, filePath):
        self.blosum.read_file(filePath)


    def distance_alignment(self, seqA, seqB):
        self.func = min
        self.scorer = self.blosum.get_distance_score
        self.gapInitCost = -self.blosum.gapInitCost
        self.gapExtendCost = -self.blosum.gapExtendCost
        self.align(seqA, seqB)

    def similarity_alignment(self, seqA, seqB):
        self.func = max
        self.scorer = self.blosum.get_similarity_score
        self.gapInitCost = self.blosum.gapInitCost
        self.gapExtendCost = self.blosum.gapExtendCost
        self.align(seqA, seqB)

    def align(self, seqA, seqB):

        self.seqA = seqA
        self.seqB = seqB

        if len(seqA) == 0 and len(seqB) == 0:
            return ()
        elif len(seqA) == 0:
            return ()
        elif len(seqB) == 0:
            return ()

        self.init_matrices()
        self.compute_matrices()

        return tuple([True, True])


    def init_matrices(self):

        self.matrixWidth = len(self.seqA) + 1
        self.matrixHeight = len(self.seqB) + 1

        self.init_alignment_matrix()
        self.init_insertion_matrix()
        self.init_deletion_matrix()


    def init_alignment_matrix(self):

        self.alignmentMatrix = [[None for j in range(self.matrixWidth)] for i in range(self.matrixHeight)]
        self.alignmentMatrix[0][0] = 0

        for i in range(1, self.matrixHeight):
            self.alignmentMatrix[i][0] = self.gapInitCost + (i * self.gapExtendCost)

        for j in range(1, self.matrixWidth):
            self.alignmentMatrix[0][j] = self.gapInitCost + (j * self.gapExtendCost)

    def init_insertion_matrix(self):

        self.insertionMatrix = [[None for j in range(self.matrixWidth)] for i in range(self.matrixHeight)]

        for i in range(self.matrixHeight):
            self.insertionMatrix[i][0] = self.alignmentMatrix[i][0] + self.gapInitCost


    def init_deletion_matrix(self):

        self.deletionMatrix = [[None for j in range(self.matrixWidth)] for i in range(self.matrixHeight)]

        for j in range(1, self.matrixWidth):
            self.deletionMatrix[0][j] = self.alignmentMatrix[0][j] + self.gapInitCost


    def compute_matrices(self):

        for i in range(1, self.matrixHeight):
            for j in range(1, self.matrixWidth):
                self.deletionMatrix[i][j] = self.func(self.deletionMatrix[i - 1][j] + self.gapExtendCost,
                                                      self.alignmentMatrix[i - 1][j] + self.gapInitCost + self.gapExtendCost)
                self.insertionMatrix[i][j] = self.func(self.insertionMatrix[i][j - 1] + self.gapExtendCost,
                                                       self.alignmentMatrix[i][j - 1] + self.gapInitCost + self.gapExtendCost)

                matchScore = self.scorer(self.seqA[j - 1], self.seqB[i - 1])
                self.alignmentMatrix[i][j] = self.func(self.alignmentMatrix[i - 1][j - 1] + matchScore,
                                                       self.insertionMatrix[i][j],
                                                       self.deletionMatrix[i][j])

