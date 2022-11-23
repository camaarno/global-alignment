#!/usr/bin/env python

"""
Contains a class representins a BLOSUM matrix. After instantiation, users may provide a file
path to the class through the read_file() method.  Data from the file is stored in the object
as a hash table, allowing for quick access to scores in the matrix.
"""

"""
@Author: Cameron Arnold
@Data: November 22nd, 2022
"""

import re


class Blosum:
    """
    Wrapper class for a BLOSUM matrix. Allows users to read in a file representing
    a BLOSUM matrix, which is stored hash table for quick access
    """

    def __init__(self):

        # Fields are populated by the read_file() method
        self.type = ""          # ex. "BLOSUM62"
        self.__headers = []     # [ Header1, Header2, ... ]
        self.__matrix = {}      # { LabelA: { LabelA: Score, LabelB: Score, ...}, LabelB: {...}, ...}
        self.gapInit = 0        # Gap initiation cost
        self.gapExtend = 0      # Gap extension cost


    def __next_line(self, file):
        """
        Reads lines from the provided file until a non-empty line is found,
        then returns the line as a list of tokens (strings)
        :param file: a reference to a File object
        :return: a list of strings, representing the next non-empty line of the file
        """
        tokens = file.readline().split()
        while not tokens:
            tokens = file.readline().split()
        return tokens


    def read_file(self, path):
        """
        Parses a provided text file representing a BLOSUM matrix and gap initiation/extension
        penalties, then uses the data to populate members of this object
        :param path: a text file containing a BLOSUM matrix, along with the gap initiation
        and gap extension costs, in the format provided in the README.md file
        """
        f = open(path, 'r')

        # Remove the title line ("BLOSUM62")
        self.type = self.__next_line(f)[0]

        # Next non-empty line contains the column headers
        self.__headers = self.__next_line(f)

        # Initialize matrix
        self.__matrix = {}
        for header in self.__headers:
            self.__matrix[header] = {}

        # Read matrix from file
        tokens = self.__next_line(f)
        rowCount = len(self.__headers)
        rowNum = 0

        while rowNum < rowCount:

            # The first string in the line is always the row label
            rowLabel = tokens.pop(0)
            if rowLabel != self.__headers[rowNum]:
                raise Exception("BLOSUM column headers must match row headers.")

            for i in range(len(tokens)):
                colLabel = self.__headers[i]
                self.__matrix[rowLabel][colLabel] = tokens[i]

                if rowLabel != colLabel:
                    self.__matrix[colLabel][rowLabel] = tokens[i]

            # Skip empty lines and increment rowNum
            tokens = self.__next_line(f)
            rowNum += 1

        # Assume the next non-empty line was "Gap_initiation = {int}"
        gapRegex = re.compile(r"^Gap_initiation = -?\d+$")
        if not gapRegex.search(' '.join(tokens)):
            raise Exception("Invalid file format.  Next line must be 'Gap_initation = {int}'")
        self.gapInit = int(tokens[2])

        # Assume the next non-empty line is "Gap_extension = {int}"
        tokens = self.__next_line(f)
        gapRegex = re.compile(r"^Gap_extension = -?\d+$")
        if not gapRegex.search(' '.join(tokens)):
            raise Exception("Invalid file format.  Next line must be 'Gap_extension = {int}'")
        self.gapExtend = int(tokens[2])

        # Finally, close the file
        f.close()


    def get_score(self, a, b):
        """
        Returns the score of amino acids a and b, as specified by the BLOSUM matrix
        :param a: a character representing an amino acid in the BLOSUM matrix
        :param b: a character representing an amino acid in the BLOSUM matrix
        :return:
        """
        return self.__matrix[a][b]

