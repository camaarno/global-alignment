


class SequenceReader:


    def read_file(self, path):
        with open(path, 'r') as f:

            seqA = f.readline().strip()
            while not seqA:
                seqA = f.readline().strip()

            seqB = f.readline().strip()
            while not seqB:
                seqB = f.readline().strip()

        return seqA, seqB
