#primer
import os, glob
import numpy as np
import pandas as pd
import Bio
from Bio.Seq import MutableSeq, Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
from typing import Tuple

def degenerate_primer(primer:'Bio.Seq.MutableSeq') -> str:
    forward = (str(primer),)
    for i in range(0, len(primer)):
        for j in range(0, len(forward)):
            primer = MutableSeq(forward[j])
            if (primer[i] == 'A') or (primer[i] == 'C') or (primer[i] == 'G') or (primer[i] == 'T'):
                pass
            else:
                forward = degenerate_primer_list(forward, primer, i, primer[i])
    return forward

def degenerate_primer_list(forward:'str', primer:'Bio.Seq.MutableSeq', i:'int', letter:'str') -> str:
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

def forward_primer_search(species:'str', forward_primer:'tuple') -> Tuple[str, str, str]:
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
        

def reverse_primer_search(species:'str', reverse_primer_set:'tuple') -> Tuple[str, str, str]:
    primer_match_query = []
    rev_primer_set = []
    
    for i in range(0,len(reverse_primer_set)):
        reverse_primer = Seq(reverse_primer_set[i])
        reverse_primer_complement = str(reverse_primer.reverse_complement())
        primer_match_query.append(species.find(reverse_primer_complement))
        rev_primer_set.append(reverse_primer_complement)
        
    if all(item == -1 for item in primer_match_query):
        return str(''), str('N/a'), str('N/a')
    else:
        for j in range(0,len(primer_match_query)):
            if primer_match_query[j] != -1:
                amplicon_segment = species[0:primer_match_query[j]+len(reverse_primer_complement)]
                rev_primer_used = rev_primer_set[j]
                reverse_primer_position = len(amplicon_segment)-len(reverse_primer_complement)
            else:
                pass
        return amplicon_segment, rev_primer_used, reverse_primer_position

def create_PCR_amplicon(core_data:'pd.DataFrame', rev_tup:'tuple', fwd_tup:'tuple') -> pd.DataFrame:
    add_on_data = []
    all_sequnces = []
    
    for item in core_data['Record id']:
        [item_rev, rev_primer_used, reverse_primer_position] = reverse_primer_search(core_data.loc[(core_data['Record id'] == item)]['16S Sequence'].item(), rev_tup)
        
        [item_amplicon, fwd_primer_used, forward_primer_position] = forward_primer_search(item_rev, fwd_tup)
        
        add_on_data.append([core_data.loc[(core_data['Record id'] == item)]['Species'].item(),
                            item,
                            fwd_primer_used, 
                            forward_primer_position, 
                            rev_primer_used, reverse_primer_position, 
                            round(GC(item_amplicon), 1), 
                            len(item_amplicon), 
                            item_amplicon])

    columns = ['Species', 'Record id', 'Forward Primer', 'forward_primer_position', 'Reverse Primer', 'reverse_primer_position', 'GC Content', 'Length of Amplicon', 'Amplicon',]
    calculated_data = pd.DataFrame(add_on_data, columns=columns)
    return calculated_data
