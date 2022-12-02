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

    NAME_DELIMITER = '-'

    def __init__(self):

        # For each list, element i represents the same individual

        self.ids = []               # list of ids
        self.proteinNames = []      # list of species names
        self.speciesNames = []      # list of species names
        self.sequences = []         # list of sequences
        self.count = 0              # total number of registered sequences


    def set_sequences(self, path):
        """
        Reads the provided text file containing sequence information, then uses the
        data to populate members of this object. Review the user manual for more
        information on the valid file format
        :param path: the file path
        :type path: str
        """
        with open(path, 'r') as f:

            # sequences come in groups of 3 lines (one line for the ID,
            # another for the name/species, and a third for the sequence)
            # keep track of which line we should be reading
            lineNum = 0

            for line in f:

                # remove unnecessary whitespaces
                line = line.strip()

                # skip empty lines
                if not line:
                    continue

                n = lineNum % 3
                if n == 0:          # First Line: ID
                    self.ids.append(line)
                elif n == 1:
                    protein, species = line.rsplit(self.NAME_DELIMITER, 1)
                    self.proteinNames.append(protein.rstrip())
                    self.speciesNames.append(species.lstrip())
                else: # n == 2
                    assert line.endswith("*")
                    self.sequences.append(line[:-1])
                    self.count += 1

                lineNum += 1

