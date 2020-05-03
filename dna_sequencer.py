import argparse
import os.path
import logging
import naive_sequencer as naive
import sequence_file_reader as reader 
import bm_preproc as bm
import pigeon_hole as ph
import overlap_finder as of
import time as t


logger = logging.getLogger()
logger.setLevel(logging.INFO)


def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--ref', action='store', default='ref.fasta',
                    help='Filename of the reference genome')

    parser.add_argument('-s', '--samp', action='store', default='sample.fasta',
                    help='Filename of the sample genome')

    parser.add_argument('-m', '--mismatchesAllowed', action='store', type=int, 
                    default=0, help='Maximum number of mismatches allowed')

    parser.add_argument('-k', '--kmerSize', action='store', type=int, 
                    default=0, help='k-mer size for pigeon hole algorithm')
    
    parser.add_argument('-i', '--ival', action='store', type=int, 
                    default=0, help='ival for subsequence search')
    
    parser.add_argument('-p', '--pigeonHole', action='store_true', 
                    help='Use pigeon hole algorithm for inexact matches')

    parser.add_argument('-b', '--boyerMoore', action='store_true', 
                    help='Use Boyer Moore algorithm for matching')

    parser.add_argument('-e', '--editDistance', action='store_true', 
                    help='Find the edit distance between the text and pattern')
    
    parser.add_argument('-c', '--reverseComplement', action='store_true', 
                    help='Also account for reverse complement with naive matching')

    parser.add_argument('-t', '--textAlphabet', action='store_true', 
                    help='Use lower-case alphabet set for Boyer Moore (default is only "ACGT"')

    parser.add_argument('-o', '--overlapReadFinder', action='store_true', 
                    help='Find overlapping reads')

    args = parser.parse_args()

    if (not os.path.isfile(args.ref)):
        print("Reference File %s does not exist" %args.ref )
        exit(-1)

    if (not os.path.isfile(args.samp)):
        print("Sample File %s does not exist" %args.samp )
        exit(-1)

    if (args.mismatchesAllowed > 0 and args.boyerMoore):
        print("Boyer Moore does not allow for mismatches. It's an exact matching algorithm")
        exit(-1)

    if (args.reverseComplement > 0 and args.boyerMoore):
        print("Reverse complement (-c) only allowed with Naive Matching (-b 0)")
        exit(-1) 

    if (args.textAlphabet > 0 and not args.boyerMoore):
        print("Alphabet set (-t) only required for Boyer Moore")
        exit(-1) 
    
    if (args.kmerSize <= 0 and args.pigeonHole):
        print("Pigeon Hole based comparisons require a k-mer size of > 0")
        exit(-1) 
    
    if (args.kmerSize <= 0 and args.overlapReadFinder):
        print("Overlap reads finder requries kmer size of > 0")
        exit(-1)

    return args


def main():
    args = parseArgs()
    ref = reader.readGenome(args.ref)
    samp = reader.readGenome(args.samp)
    mismatchesAllowed = args.mismatchesAllowed
    logging.debug("Sample: " + samp)
    logging.debug("Reference: " + ref)
    matches = []
    alignments  = 0
    numComparisons = 0
    queries = 0
    editDistance = 0
    overlaps = []
    duration = 0
    if (not (args.boyerMoore or args.pigeonHole or args.editDistance or args.overlapReadFinder)):
        # Naive matching algorithm
        if (args.reverseComplement):
            matches = naive.naiveWithReverseComplement(samp, ref, mismatchesAllowed)
        else:
            start = t.perf_counter()
            matches, alignments, numComparisons = naive.naiveWithNMismatches(samp, ref, mismatchesAllowed)
            duration = t.perf_counter() - start
    else:
        if ( args.boyerMoore):
            if (not args.textAlphabet):
                p_bm = bm.BoyerMoore(samp, bm.DNA_SEQUENCING_ALPHABET)
            else: 
                p_bm = bm.BoyerMoore(samp, bm.LOWERCASE_ALPHABET)
            matches, alignments, numComparisons = bm.boyer_moore(samp, p_bm, ref)

        if (args.pigeonHole):
            # Pigeon Hole 
            # t, p, k, m, i=0):
            p_ph = ph.PigeonHole(ref, samp, args.kmerSize, args.mismatchesAllowed, args.ival)
            start = t.perf_counter()
            matches, queries, numComparisons = p_ph.getMatches()
            duration = t.perf_counter() - start

        if (args.editDistance):
            start = t.perf_counter()
            editDistance = of.approximateMatch(ref, samp)
            duration = t.perf_counter() - start

        if (args.overlapReadFinder):
            reads, _ = reader.readFastq(args.samp)
            start = t.perf_counter()
            overlaps = of.findOverlapsForReads(reads, args.kmerSize)
            duration = t.perf_counter() - start
            sources = of.findSourceNodes(overlaps)
            logging.debug("All overlaps: %s" %overlaps)

    print ("Matches                           :   %s, Total matches: %d" %(matches, len(matches)))
    print ("Aligments (naive and Boyer Moore) :   %d" %(alignments))
    print ("Queries (pigeon-hole only)        :   %d" %(queries))
    print ("Char comparisons                  :   %d" %(numComparisons))
    print ("Duration (seconds)                :   %f" %(duration))
    print ("Edit Distance (edit-distance only):   %d" %(editDistance))
    print ("Num Overlaps (overlap finder only):   %s" %(len(overlaps)))
    print ("Num sources (overlap finder only) :   %s" %(len(sources)))


if __name__ == "__main__":
    main()
