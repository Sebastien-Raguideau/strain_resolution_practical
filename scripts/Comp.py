
from itertools import compress
import argparse
import sys
import numpy as np
import os
import subprocess
import re 

from numpy.random import RandomState
import logging 

from subprocess import PIPE
from collections import defaultdict
from itertools import product, tee
from collections import Counter, OrderedDict
from sklearn.decomposition import PCA
from sklearn import preprocessing
import scipy as sp

def readFasta(fastaFileName):

    fasta = defaultdict(str)
    ids = []
    lens = {}
    
    with open(fastaFileName) as file_one:
        for line in map(str.rstrip, file_one):
            if line.startswith(">"):
                sequence_name = line.lstrip(">")
                ids.append(sequence_name)
            else:
                fasta[sequence_name] += line
    
    for iid in ids:
        lens[iid] = len(fasta[iid])
    
    return (ids,fasta,lens)

def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in range(i):
            next(el, None)
    return zip(*els)

def _calculate_composition(fasta, length_threshold, kmer_len):
    #Generate kmer dictionary
    feature_mapping, nr_features = generate_feature_mapping(kmer_len)

    # Store composition vectors in a dictionary before creating dataframe
    composition_d = OrderedDict()
    contig_lengths = OrderedDict()
    
    for (sid,seq) in fasta.items():
    
        seq_len = len(seq)
    
        if seq_len<= length_threshold:
            continue
            
        contig_lengths[sid] = seq_len
        # Create a list containing all kmers, translated to integers
        kmers = [
                feature_mapping[kmer_tuple]
                for kmer_tuple 
                in window(seq.upper(), kmer_len)
                if kmer_tuple in feature_mapping
                ]
        # numpy.bincount returns an array of size = max + 1
        # so we add the max value and remove it afterwards
        # numpy.bincount was found to be much more efficient than
        # counting manually or using collections.Counter
        kmers.append(nr_features - 1)
        composition_v = np.bincount(np.array(kmers))
        composition_v[-1] -= 1
        # Adding pseudo counts before storing in dict
        
        composition_d[sid] = composition_v + np.ones(nr_features)
    

    return composition_d, contig_lengths, nr_features



def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC",repeat=kmer_len):
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = tuple([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash, counter
    
def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)
        
        


        
def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument("fastaFile", help="gfa file")
    
    args = parser.parse_args()

   # import ipdb; ipdb.set_trace()    

    (ids,fasta,lens) = readFasta(args.fastaFile)
    
    (composition_d, contig_lengths,nr_features) = _calculate_composition(fasta, 10, 4)
    
    N = len(composition_d)
    
    comp_array = np.zeros((N,nr_features))
    
    uMap = {}
    sids = []
    sidx = 0
    for (sid,comp) in composition_d.items():
    
        comp_array[sidx,:] = comp
        
        sids.append(sid)
        
        uMap[sid] = sidx
        
        sidx += 1
    
    
    row_sums = comp_array.sum(axis=1)
    comp_arrayP =  np.log(comp_array / row_sums[:, np.newaxis])
    
    pca_object = PCA()
    
    pca_object.fit(comp_arrayP)
    
    tran_comp = pca_object.transform(comp_arrayP)
    
    
    D = np.argmax(np.cumsum(pca_object.explained_variance_ratio_) > 0.9)
    
    with open('PCA_tran2.csv','w') as f:
        
        hString = ','.join([ 'D' + str(x) for x in range(D)])
        print('Unitig,%s' % (hString),file=f)
        
        for u, uid in enumerate(sids):
            xString = ','.join([str(x) for x in (tran_comp[u,0:D]).tolist()])
            
            print('%s,%s' % (uid,xString),file=f)
            

    
    
    
    
    min_max_scaler = preprocessing.MinMaxScaler()
    tran_comp_minmax = min_max_scaler.fit_transform(tran_comp)
    print('NODE,COLOUR') 
    for (s,sid) in enumerate(sids):
    
        r = int(255*tran_comp_minmax[s,0])
        if r == 255:
            r -= 1
        
        g = int(255*tran_comp_minmax[s,1])
        if g == 255:
            g -= 1
            
        b = int(255*tran_comp_minmax[s,2])
        if b == 255:
            b -= 1
             
        hexCol = rgb_to_hex(r,g,b)
        
        print('%s,%s' % (sid,hexCol))
    

    
if __name__ == "__main__":
    main(sys.argv[1:])
    
    
