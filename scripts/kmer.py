# from gotoh2 import *
import argparse
import csv
import math
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate k-mer distance matrix for all protein sequences "
                    "in a FASTA file."
    )
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='input, FASTA file with protein sequences')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help='output, distance matrix as CSV with accession and '
                             'gene name as column names.')
    parser.add_argument('-header', type=argparse.FileType('w'), default=None,
                        help='output (optional), protein header info as CSV')
    parser.add_argument('--kernel', action='store_true',
                        help='optional, use string kernel (d2s) instead of '
                        'default intersection distance.')
    parser.add_argument('--filter', action='store_true',
                        help='<optional> Get a list with accession numbers'
                         'and number of proteins associated with each accession.')
    return parser.parse_args()


def kmer(seq):
    """
    Calculate word counts for word lengths varying from k={1,2,3}
    :param seq:  str, protein sequence
    :return:  dict, counts keyed by word - absent entries imply zero
    """
    d = {}
    for i in range(len(seq)):
        for k in [1, 2, 3]:
            if (len(seq)-1 - i) < k:
                break
            word = seq[i:(i+k)]
            if word not in d:
                d.update({word: 0})
            d[word] += 1
    return d


def kdist(k1, k2):
    """
    string kernel / d2s distance
    :param k1:  kmer counts for sequence 1
    :param k2:  kmer counts for sequence 2
    :return:  pairwise distance
    """
    d11, d12, d22 = 0, 0, 0
    for km in set(k1.keys()).union(set(k2.keys())):
        c1 = k1.get(km, 0)
        c2 = k2.get(km, 0)
        d11 += c1 * c1
        d12 += c1 * c2
        d22 += c2 * c2

    return d12 / math.sqrt(d11 * d22)


def intersection(k1, k2):
    """
    Intersection distance
    :param k1:  kmer counts for sequence 1
    :param k2:  kmer counts for sequence 2
    :return:  float, pairwise distance
    """
    res = 0
    for km in set(k1.keys()).intersection(set(k2.keys())):
        res += 2 * min(k1[km], k2[km])
    return res / (sum(k1.values()) + sum(k2.values()))


def iter_fasta(handle):
    """
    Function from gotoh2: https://github.com/ArtPoon/gotoh2/blob/master/gotoh2_utils.py
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.

    :param handle:  open stream to FASTA file in read mode
    :yield tuples, (header, sequence)
    """
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                yield h, sequence
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    yield h, sequence

def get_accns(labels):
    """
    Count number of proteins per accesssion number
    """
    accns = {}
    for label in labels:
        acc = label.split('-')[0]
        prot = label.split('-')[2]
        if acc not in accns.keys():
            accns[acc] = {'names': [], 'count':0, 'labels':[]}
        accns[acc]['names'].append(prot)
        accns[acc]['labels'].append(label)
        accns[acc]['count'] += 1
    
    return accns

def parse_fasta(infile, filter):
    """
    Parse fasta file to get sequence and headers
    If filter, remove sequences with less proteins than expected
    """
    # parse FASTA input
    acc_counts = {}
    for h, s in iter_fasta(infile):
        acc = h.split('-')[0]
        prot = h.split('-')[2]
        
        if acc not in acc_counts.keys():
            acc_counts[acc] = {'names': [], 'count':0, 
                               'label':[], 'seq':[]}
        acc_counts[acc]['names'].append(prot)
        acc_counts[acc]['label'].append(h)
        acc_counts[acc]['seq'].append(s)
        acc_counts[acc]['count'] += 1
    
    # Filter accession numbers below a given number of proteins
    labels = []
    seqs = []
    counts = [acc_counts[accn]['count'] for accn in acc_counts.keys()]
    mean = sum(counts) / len(counts)
    print(f'Mean number of proteins per genome: {mean:1.2f}\n')
    ignored = 0
    for accn in acc_counts.keys():
        n_prots = acc_counts[accn]['count']
        if filter and n_prots < (mean-2):
                print (f'Skip: {accn}. {n_prots} proteins: {acc_counts[accn]["names"]}')
                ignored += 1
                continue
                # print(f"Protein names: {acc_counts[accn]['names']}")     
        labels.extend(acc_counts[accn]['label'])
        seqs.extend(acc_counts[accn]['seq'])
    print(f'accession numbers ignored: {ignored}')
    return labels, seqs

def kmer_dist(infile, outfile, header, kernel, filter):
    """
    :param infile:  open file stream to read FASTA
    :param outfile:  open file stream to write distances as CSV
    :param header:  open file stream to write header info as CSV
    :param kernel:  if True, calculate kernel distance instead of
                    intersection distance
    """

    labels, seqs = parse_fasta(infile, filter)

    # pre-calculate k-mer counts
    kmers = {}
    for i, seq in enumerate(seqs):
        kmers.update({labels[i]: kmer(seq)})

    # calculate and write distance matrix
    outstr = None
    for i, l1 in enumerate(labels):
        if i > 0:
            sys.stdout.write('\b' * len(outstr))

        outstr = f'{i} / {len(labels)}'
        sys.stdout.write(outstr)
        sys.stdout.flush()

        outfile.write(f'"{l1}"')
        for l2 in labels:
            if kernel:
                d = 1. if l1 == l2 else kdist(kmers[l1], kmers[l2])
            else:
                # default
                d = intersection(kmers[l1], kmers[l2])
            
            outfile.write(f',{d:1.9f}')

        outfile.write('\n')

    sys.stdout.write('\n')


if __name__ == '__main__':
    args = parse_args()
    kmer_dist(infile=args.infile, outfile=args.outfile, header=args.header,
             kernel=args.kernel, filter = args.filter)
