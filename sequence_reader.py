#!/usr/bin/env python

"""
Contains the SequenceReader class
"""

"""
@Author: Cameron Arnold
@Data: November 27th, 2022
"""


class SequenceReader:
    """
    A class for reading sequence files, whereas a sequence file is a .txt file
    with the first two non-empty lines each containing an amino acid sequence.
    """

    def read_file(self, path):
        """
        Reads and returns the two sequences in the sequence file at path
        :param path: the file path
        :type path: str
        :return: a tuple containing the two sequences
        :rtype: tuple
        """
        with open(path, 'r') as f:

            seqA = f.readline().strip()
            while not seqA:   # ignore empty lines
                seqA = f.readline().strip()

            seqB = f.readline().strip()
            while not seqB:   # ignore empty lines
                seqB = f.readline().strip()

        return seqA, seqB
