import kmer_index as idx

class PigeonHole(object):
    """ Uses the Index object to find all strings with up to 'n' mismatches"""

    def __init__(self, t, p, k, m, i=0):
        """ Create an Index of kmer-length k. Store the pattern and string
            t: text
            p: pattern
            k: k-mer
            m: mismatches allowed
            i: ival -- use subsequence search
        """
        self.t = t
        self.p = p
        self.k = k
        self.m = m
        self.ival = i
        if (i == 0):
            self.idx = idx.Index(t,k)
        else:
            self.idx = idx.SubseqIndex(t, k, i)
        
    def naiveWithNMisMatches(self, p, t):
        numMismatches = 0
        numComparisons = 0
        for j in range(len(p)):  # loop over characters
            numComparisons += 1
            if t[j] != p[j]:  # compare characters
                numMismatches += 1
                if numMismatches > self.m:
                    return False, numComparisons
        return True, numComparisons

    def getMatches(self):
        if (self.ival == 0):
            numSlices = int(len(self.p) / self.k) + 1
        else:
            numSlices = self.ival 
        matches = set()
        queries = 0
        numComparisons = 0
        for s in range(numSlices):
            hits = self.idx.query(self.p, s)
            for hit in hits:
                queries += 1
                if (self.ival == 0):
                    start = hit - s*self.k
                else:
                    start = hit - s
                if ( not start in matches):
                    match, comparisons = self.naiveWithNMisMatches(self.p, self.t[start: start + len(self.p)])
                    if (match):
                        matches.add(start)
                    numComparisons += comparisons
        return sorted(matches), queries, numComparisons
        