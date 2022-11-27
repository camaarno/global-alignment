#!/usr/bin/env python

"""
Contains the GlobalAligner class
"""

"""
@Author: Cameron Arnold
@Data: November 27th, 2022
"""


from blosum_reader import BlosumReader  # used to read BLOSUM file
from operator import lt, gt             # used to compare scores


class GlobalAligner:
    """
    Computes similarity and distance alignments on provided sequences, using
    the affine indel gap model fit sequence alignment
    """

    ALIGNMENT_SYMBOL = "A"
    DELETION_SYMBOL = "D"
    INSERTION_SYMBOL = "I"
    GAP_SYMBOL = '-'

    def __init__(self):

        self.seqX = ""                  # horizontal sequence
        self.seqY = ""                  # vertical sequence

        self.blosum = BlosumReader()    # BLOSUM file reader
        self.gapInitCost = 0            # gap initiation cost
        self.gapExtendCost = 0          # gap extension cost

        self.alignmentMatrix = []       # optimal alignment matrix (A)
        self.deletionMatrix = []        # alignment ends with deletion (D)
        self.insertionMatrix = []       # alignment ends with insertion (I)

        self.matrixHeight = 0           # height of each matrix
        self.matrixWidth = 0            # width of each matrix

        self.compare_function = None    # lt / gt
        self.score_function = None      # blosum.get_similarity_score / blosum.get_distance_score


    def set_blosum(self, path):
        """
        Get the BLOSUM matrix and gap costs from the file at the provided path.
        See user manual for more information on the file format required.
        :param path: the file path
        :type path: str
        """
        self.blosum.read_file(path)


    def distance_alignment(self, seqX, seqY):
        """
        Compute the global distance alignment of seqX and seqY
        :param seqX: an amino acid sequence
        :type seqX: str
        :param seqY: an amino acid sequence
        :type seqY: str
        :return: an aligned seqX, an aligned seqY, and the score
        :rtype: tuple
        """
        self.compare_function = lt
        self.score_function = self.blosum.get_distance_score
        self.gapInitCost = -self.blosum.gapInitCost
        self.gapExtendCost = -self.blosum.gapExtendCost
        return self.align(seqX, seqY)

    def similarity_alignment(self, seqX, seqY):
        """
        Compute the global simlarity alignment of seqX and seqY
        :param seqX: an amino acid sequence
        :type seqX: str
        :param seqY: an amino acid sequence
        :type seqY: str
        :return: an aligned seqX, an aligned seqY, and the score
        :rtype: tuple
        """
        self.compare_function = gt
        self.score_function = self.blosum.get_similarity_score
        self.gapInitCost = self.blosum.gapInitCost
        self.gapExtendCost = self.blosum.gapExtendCost
        return self.align(seqX, seqY)

    def align(self, seqX, seqY):
        """
        Compute the alignment of seqX and seqY.  Requires that compare_function,
        score_function, gapInitCost, and gapExtend cost are set to a valid value
        :param seqX: an amino acid sequence
        :type seqX: str
        :param seqY: an amino acid sequence
        :type seqY: str
        :return: an aligned seqX, an aligned seqY, and the score
        :rtype: tuple
        """

        if not seqX or not seqY:
            return seqX, seqY

        self.seqX = seqX
        self.seqY = seqY

        self.init_matrices()
        self.compute_matrices()
        result = self.traceback_matrices()
        return result


    def init_matrices(self):
        """
        Initialize the three matrices
        """
        self.matrixWidth = len(self.seqX) + 1
        self.matrixHeight = len(self.seqY) + 1

        self.init_alignment_matrix()
        self.init_insertion_matrix()
        self.init_deletion_matrix()


    def init_alignment_matrix(self):
        """
        Initialize the alignment matrix (A) by allocating space for the matrix and setting
        the values for the first row and the first column
        """

        self.alignmentMatrix = [[None for j in range(self.matrixWidth)] for i in range(self.matrixHeight)]
        self.alignmentMatrix[0][0] = (0, None)

        for i in range(1, self.matrixHeight):
            self.alignmentMatrix[i][0] = (self.gapInitCost + (i * self.gapExtendCost), None)

        for j in range(1, self.matrixWidth):
            self.alignmentMatrix[0][j] = (self.gapInitCost + (j * self.gapExtendCost), None)


    def init_insertion_matrix(self):
        """
        Initialize the insertion matrix (I) by allocating space for the matrix, and setting
        the values of the first column
        """

        self.insertionMatrix = [[None for j in range(self.matrixWidth)] for i in range(self.matrixHeight)]

        for i in range(self.matrixHeight):
            self.insertionMatrix[i][0] = (self.alignmentMatrix[i][0][0] + self.gapInitCost, None)


    def init_deletion_matrix(self):
        """
        Initialize the deletion matrix (D) by allocating space for the matrix, and setting
        the values of the first row
        """

        self.deletionMatrix = [[None for j in range(self.matrixWidth)] for i in range(self.matrixHeight)]

        for j in range(1, self.matrixWidth):
            self.deletionMatrix[0][j] = (self.alignmentMatrix[0][j][0] + self.gapInitCost, None)


    def compute_matrices(self):
        """
        Populate the three matrices using the affine indel gap model global sequence
        alignment algorithm
        """

        #costD = None
        #costI = None
        #costA = None

        for i in range(1, self.matrixHeight):
            for j in range(1, self.matrixWidth):

                costD = self.deletionMatrix[i - 1][j][0] + self.gapExtendCost
                costA = self.alignmentMatrix[i - 1][j][0] + self.gapInitCost + self.gapExtendCost
                self.deletionMatrix[i][j] = (costD, None) if self.compare_function(costD, costA) else (costA, self.ALIGNMENT_SYMBOL)

                costI = self.insertionMatrix[i][j - 1][0] + self.gapExtendCost
                costA = self.alignmentMatrix[i][j - 1][0] + self.gapInitCost + self.gapExtendCost
                self.insertionMatrix[i][j] = (costI, None) if self.compare_function(costI, costA) else (costA, self.ALIGNMENT_SYMBOL)

                matchScore = self.score_function(self.seqX[j - 1], self.seqY[i - 1])
                costA = self.alignmentMatrix[i - 1][j - 1][0] + matchScore
                costI = self.insertionMatrix[i][j][0]
                costD = self.deletionMatrix[i][j][0]

                if self.compare_function(costA, costI):
                    self.alignmentMatrix[i][j] = (costA, None) if self.compare_function(costA, costD) else (costD, self.DELETION_SYMBOL)
                else:
                    self.alignmentMatrix[i][j] = (costI, self.INSERTION_SYMBOL) if self.compare_function(costI, costD) else (costD, self.DELETION_SYMBOL)

    def traceback_matrices(self):
        """
        Traceback the three matrices in order to compute the alignment
        :return: an aligned seqX, an aligned seqY, and the score
        :rtype: tuple
        """

        i = self.matrixHeight - 1
        j = self.matrixWidth - 1
        currentMatrix = self.ALIGNMENT_SYMBOL
        matrix = self.alignmentMatrix

        alignX = []
        alignY = []

        while i > 0 and j > 0:

            if currentMatrix == self.ALIGNMENT_SYMBOL:
                matrix = self.alignmentMatrix
            elif currentMatrix == self.INSERTION_SYMBOL:
                matrix = self.insertionMatrix
            elif currentMatrix == self.DELETION_SYMBOL:
                matrix = self.deletionMatrix

            _, originMatrix = matrix[i][j]

            if currentMatrix == self.ALIGNMENT_SYMBOL:
                if originMatrix is None:  # same matrix (A)
                    alignX.append(self.seqX[j - 1])
                    alignY.append(self.seqY[i - 1])
                    j -= 1
                    i -= 1
                elif originMatrix == self.INSERTION_SYMBOL:
                    currentMatrix = self.INSERTION_SYMBOL
                elif originMatrix == self.DELETION_SYMBOL:
                    currentMatrix = self.DELETION_SYMBOL

            elif currentMatrix == self.INSERTION_SYMBOL:
                if originMatrix is self.ALIGNMENT_SYMBOL:
                    currentMatrix = self.ALIGNMENT_SYMBOL
                alignX.append(self.seqX[j - 1])
                alignY.append(self.GAP_SYMBOL)
                j -= 1

            elif currentMatrix == self.DELETION_SYMBOL:
                if originMatrix is self.ALIGNMENT_SYMBOL:
                    currentMatrix = self.ALIGNMENT_SYMBOL
                alignX.append(self.GAP_SYMBOL)
                alignY.append(self.seqY[i - 1])
                i -= 1

        alignX = ''.join(alignX[::-1])
        alignY = ''.join(alignY[::-1])
        score = self.alignmentMatrix[self.matrixHeight - 1][self.matrixWidth - 1][0]

        return alignX, alignY, score

