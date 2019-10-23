import pandas as pd
from Bio import SeqIO
import argparse
from multiprocessing import Pool
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fastafile'    , help='file to process, contigs fasta'          , type=str, action = 'store', required = True)
parser.add_argument('-k', '--kmersize'     , help='Kmer size to use for analysis'           , type=int, action = 'store', required = True)
parser.add_argument('-t', '--ThreadNumber' , help='Number of thread for parrallel computing', type=int, action = 'store', required = True)

args = parser.parse_args()
fasta_file = args.fastafile
k = args.kmersize
t = args.ThreadNumber

# Get Kmer spectrum of a sequence
def KmersFreqSeq(sequence, k):
    kmer_dict = dict()
    for char in range(0,len(sequence) - k + 1):
        kmer = sequence[char:char + k]
        if 'N' in kmer:
            continue
        else:
            if kmer in kmer_dict.keys():
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    count = sum([kmer_dict[i] for i in kmer_dict.keys()])
    for kmer in kmer_dict.keys():
        kmer_dict[kmer] = kmer_dict[kmer] / count
    return(kmer_dict)

def GetGC(sequence):
    return((sequence.count('G') + sequence.count('C')) / (sequence.count('A') + sequence.count('C') + sequence.count('G') + sequence.count('T')))

def Process_file(record):
    dataset = pd.DataFrame()
    for sequence in record:
        # Calculate tnfs
        record_output = dict()
        record_output = KmersFreqSeq(str(sequence.id), k)
        # Calculate sequence length and gc
        record_output['contig'] = str(sequence.id)
        record_output['length'] = len(str(sequence.seq))
        record_output['gc_content'] = GetGC(str(sequence.seq))
        dataset = dataset.append(record_output, ignore_index=True)
    return(dataset)

def Main():
    print('Parsing fasta file...')
    main_table = pd.DataFrame()
    data = list(SeqIO.parse(fasta_file, 'fasta'))
    record_list = np.array_split(data,t)

    print('Calculating TNFs...')
    p = Pool(t)
    output = p.map(Process_file, record_list)

    print('Merging results...')
    dataset = pd.concat(output)
    dataset.to_csv(fasta_file + '_' + str(k) + 'mers.csv', sep='\t')

    print('Done!')

if __name__ == '__main__':
    Main()
