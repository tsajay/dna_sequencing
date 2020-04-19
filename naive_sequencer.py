import argparse
import os.path
import logging

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def naiveWithNMismatches(p, t, n):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        numMismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                numMismatches += 1
                if numMismatches > n:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def naiveWithReverseComplement(p, t):
    r = reverseComplement(p)
    logging.debug("rc: " + r)
    occurences = naive(p, t)
    if (p != r):
        occurences += naive(r, t)
    occurences.sort()
    return occurences
    

def naiveWithReverseComplement2(p, t):
    r = reverseComplement(p)
    logging.debug("rc: " + r)
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
        match = True
        if not p == r:
            for j in range(len(p)):  # loop over characters
                if t[i+j] != r[j]:  # compare characters
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--ref', action='store', required='true',
                    help='Filename of the reference genome')

    parser.add_argument('-s', '--samp', action='store', required='true',
                    help='Filename of the sample genome')

    parser.add_argument('-m', '--mismatchesAllowed', action='store', type=int, 
                    default=0, help='Maximum number of mismatches allowed')

    args = parser.parse_args()

    if (not os.path.isfile(args.ref)):
        print("Reference File %s does not exist" %args.ref )
        exit(-1)

    if (not os.path.isfile(args.samp)):
        print("Sample File %s does not exist" %args.samp )
        exit(-1)

    return args


def main():
    args = parseArgs()
    ref = readGenome(args.ref)
    samp = readGenome(args.samp)
    mismatchesAllowed = args.mismatchesAllowed
    logging.debug("Sample: " + samp)
    logging.debug("Reference: " + ref)
    if (mismatchesAllowed == 0):
        matches = naiveWithReverseComplement(samp, ref)
    else:
        matches = naiveWithNMismatches(samp, ref, mismatchesAllowed)
    print ("Matches: %s, Total matches: %d" %(matches, len(matches)))



if __name__ == "__main__":
    main()
