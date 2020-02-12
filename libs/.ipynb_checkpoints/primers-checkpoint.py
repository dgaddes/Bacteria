#primers

import os
import glob
import numpy as np
import pandas as pd
import Bio
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC

def degenerate_primer(primer):
    forward = (str(primer),)
    for i in range(0, len(primer)):
        for j in range(0, len(forward)):
            primer = MutableSeq(forward[j])
            if (primer[i] == 'A') or (primer[i] == 'C') or (primer[i] == 'G') or (primer[i] == 'T'):
                pass
            else:
                forward = convert_degenerate(forward, primer, i, primer[i])
    return forward

def convert_degenerate(forward, primer, i, letter):
    R = ['A', 'G', 'R']
    M = ['A', 'C', 'M']
    S = ['C', 'G', 'S']
    B = ['C', 'G', 'T', 'B']
    H = ['A', 'C', 'T', 'H']
    N = ['A', 'C', 'G', 'T', 'N']
    Y = ['C', 'T', 'Y']
    K = ['G', 'T', 'K']
    W = ['A', 'T', 'W']
    D = ['A', 'G', 'T', 'D']
    V = ['A', 'C', 'G', 'V']
    mixed_nucleotides = [R, M, S, B, H, N, Y, K, W, D, V]
    mixed_strings = ['R', 'M', 'S', 'B', 'H', 'N', 'Y', 'K', 'W', 'D', 'V']
    k = 0
    for string in mixed_strings:
        if letter == string:
            break
        else:
            k = k+1
            
    for basepair in mixed_nucleotides[k]:
        primer[i] = basepair
        forward = forward + (str(primer),)
    return forward

def forward_primer_search(species, forward_primer):
    primer_match_query = []
    fwd_primer_set = []
    init_len = len(species)
    for i in range(0,len(forward_primer)):
        primer_match_query.append(species.find(forward_primer[i]))
        fwd_primer_set.append(forward_primer[i])
    
    if all(item == -1 for item in primer_match_query):
        return str(''), str('N/a'), str('N/a')
    
    else:
        for k in range(0, len(primer_match_query)):
            if primer_match_query[k] != -1:
                forward_amplicon_segment = species[primer_match_query[k]:len(species)]
                fwd_primer_used = forward_primer[k]
                foward_primer_position = len(species) - len(forward_amplicon_segment)
            else:
                pass
        return forward_amplicon_segment, fwd_primer_used, foward_primer_position
        

def reverse_primer_search(species, reverse_primer):
    primer_match_query = []
    rev_primer_set = []
    for i in range(0,len(reverse_primer)):
        rev = Seq(reverse_primer[i])
        rev = str(rev.reverse_complement())
        primer_match_query.append(species.find(rev))
        rev_primer_set.append(rev)
    if all(item == -1 for item in primer_match_query):
        return str(''), str('N/a'), str('N/a')
    else:
        for j in range(0,len(primer_match_query)):
            if primer_match_query[j] != -1:
                amplicon_segment = species[0:primer_match_query[j]+len(rev)]
                rev_primer_used = rev_primer_set[j]
                reverse_primer_position = len(species) - len(amplicon_segment)
            else:
                pass
        return amplicon_segment, rev_primer_used, reverse_primer_position

def create_PCR_amplicon(core_data, rev_tup, fwd_tup):
    add_on_data = []
    all_sequnces = []
    for item in core_data['Species']:
        #a = core_data.loc[core_data['Species'] == item]
        
        [item_rev, rev_primer_used, reverse_primer_position] = reverse_primer_search(core_data.loc[(core_data['Species'] == item)]['16S Sequence'].item(), rev_tup)
        
        [item_amplicon, fwd_primer_used, forward_primer_position] = forward_primer_search(item_rev, fwd_tup)
        
        add_on_data.append([item, fwd_primer_used, forward_primer_position, rev_primer_used, reverse_primer_position, round(GC(item_amplicon), 1), len(item_amplicon), item_amplicon])

    columns = ['Species', 'Forward Primer', 'forward_primer_position', 'Reverse Primer', 'reverse_primer_position', 'GC Content', 'Length of Amplicon', 'Amplicon',]
    calculated_data = pd.DataFrame(add_on_data, columns=columns)
    return calculated_data
