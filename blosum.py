#!/usr/bin/env python

"""
Contains a class representing a BLOSUM matrix. After instantiation, users may provide a file
path to the class through the read_file() method.  Data from the file is stored in the object
as a hash table, allowing for quick access to scores in the matrix.
"""

"""
@Author: Cameron Arnold
@Data: November 23rd, 2022
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
        self._headers = []     # [ Header1, Header2, ... ]
        self._matrix = {}      # { LabelA: { LabelA: Score, LabelB: Score, ...}, LabelB: {...}, ...}
        self.gapInit = 0        # Gap initiation cost
        self.gapExtend = 0      # Gap extension cost


    def _next_line(self, file):
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
        with open(path, 'r') as f:
            self._read_headers(f)
            self._create_matrix()
            nextLine = self._read_matrix(f)
            self._read_gaps(f, nextLine)


    def _read_headers(self, f):
        """
        Reads the headers from file f, storing it as member _headers
        :param f: the file from which to read the headers
        """
        self.type = self._next_line(f)[0]       # Remove the title line (ex. "BLOSUM62")
        self._headers = self._next_line(f)      # Next non-empty line contains the column headers


    def _create_matrix(self):
        """
        Create the initial matrix, using the header labels in member _headers.
        """
        self._matrix = {}
        for header in self._headers:
            self._matrix[header] = {}


    def _read_matrix(self, f):
        """
        Reads the matrix from file f, storing the data in the hashtable member _matrix
        :param f: the file from which to read the matrix
        :return: a list of strings, representing the "next" non-empty line in the file
        """
        tokens = self._next_line(f)
        rowCount = len(self._headers)
        rowNum = 0

        while rowNum < rowCount:

            # The first string in the line is always the row label
            rowLabel = tokens.pop(0)
            if rowLabel != self._headers[rowNum]:
                raise Exception("BLOSUM column headers must match row headers.")

            for i in range(len(tokens)):
                colLabel = self._headers[i]
                self._matrix[rowLabel][colLabel] = tokens[i]

                if rowLabel != colLabel:
                    self._matrix[colLabel][rowLabel] = tokens[i]

            # Skip empty lines and increment rowNum
            tokens = self._next_line(f)
            rowNum += 1

        return tokens


    def _read_gaps(self, f, firstLine):
        """
        Reads the gap initiation and gap extension costs from file f, and stores the
        data in members gapInit and gapExtend
        :param f: the file from which to read the gap costs
        :param firstLine: a list of strings, representing the line containing the gap initiation cost
        """
        # Assume the next non-empty line was "Gap_initiation = {int}"
        gapRegex = re.compile(r"^Gap_initiation = -?\d+$")
        if not gapRegex.search(' '.join(firstLine)):
            raise Exception("Invalid file format.  Next line must be 'Gap_initation = {int}'")
        self.gapInit = int(firstLine[2])

        # Assume the next non-empty line is "Gap_extension = {int}"
        tokens = self._next_line(f)
        gapRegex = re.compile(r"^Gap_extension = -?\d+$")
        if not gapRegex.search(' '.join(tokens)):
            raise Exception("Invalid file format.  Next line must be 'Gap_extension = {int}'")
        self.gapExtend = int(tokens[2])


    def get_score(self, a, b):
        """
        Returns the score of amino acids a and b, as specified by the BLOSUM matrix
        :param a: a character representing an amino acid in the BLOSUM matrix
        :param b: a character representing an amino acid in the BLOSUM matrix
        :return: the BLOSUM score of amino acids a and b
        """
        return self._matrix[a][b]

