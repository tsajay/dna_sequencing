import kmer_index as idx
'''
import bisect


class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
'''

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
                    #if (self.naiveWithNMisMatches(self.p, self.t[start: start + len(self.p)]):
                    match, comparisons = self.naiveWithNMisMatches(self.p, self.t[start: start + len(self.p)])
                    if (match):
                        matches.add(start)
                    numComparisons += comparisons
        return sorted(matches), queries, numComparisons

                
        



        
    