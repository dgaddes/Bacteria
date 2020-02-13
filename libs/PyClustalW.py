#PyClustalW
'''
__author__ = "David Gaddes"
__Copyright__ "Copyright September 2019, David Gaddes"
__License__ = "GPL"
__email__ "dgaddes@protonmail.com"  
'''

import pandas as pd
import Bio
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo,SeqIO

def create_alignment_file(pd_frame:'pd.DataFrame', column_name:'str', file_name:'str') -> Bio.Align.MultipleSeqAlignment:
    '''This function creates a fasta file suitable for alignment'''
    
    alignment_list = []
    min_legnth = pd_frame[column_name].map(lambda x: len(x) if len(x) > 0 else 10**5).min()
    
    for item in pd_frame['Record id']:
        if len(pd_frame.loc[(pd_frame['Record id'] == item)][column_name].item()) > 1:
            alignment_list.append(SeqRecord(Seq(pd_frame.loc[(pd_frame['Record id'] == item)][column_name].item()[1:min_legnth], generic_dna), 
            id = item,
            description = pd_frame.loc[(pd_frame['Record id'] == item)]['Description'].item()))
        else:
            pass
        
    writen_alignment = MultipleSeqAlignment(alignment_list)
    AlignIO.write(writen_alignment, str(file_name +'.faa'), "fasta")
    return writen_alignment
    
def ClustalW_alignment(writen_alignment:'str') -> None:
    '''This function runs the alignment from the command line'''
    clustalw_cline = ClustalwCommandline(infile=writen_alignment)
    stdout, stderr = clustalw_cline()

def Phylo_tree(phylo_str:'str') -> None:
    '''This function generates a phylotree'''
    tree = Phylo.read(str(phylo_str+'.dnd'), 'newick')
    Phylo.draw_ascii(tree)
    
def Read_alignment_file(Read_align_str:'str') -> pd.DataFrame:
    '''This function adds the aligned data to a Pandas Dataframe'''
    alignment_record = []
    with open(str(Read_align_str+'.aln'), "rU") as handle:
        for seq_record in SeqIO.parse(handle, "clustal"):
            offset = count_alignment_offset(str(seq_record.seq))
            alignment_record.append([seq_record.id, str(seq_record.seq), offset])

    columns = ['Record id', 'Aligned 16S Sequence', 'offset']
    align_data = pd.DataFrame(alignment_record, columns=columns)
    return align_data

def count_alignment_offset(aligned_sequence:'str') -> int:
    count = 0
    if len(aligned_sequence)>1:
        for i in range(1,len(aligned_sequence)):
            if aligned_sequence[i-1]=='-':
                count+=1
            else:
                break
    return count

def conserved_regions(pd_frame:'pd.DataFrame', column:'str') -> pd.DataFrame:
    conserved_region = []
    for letter in range(0, len(pd_frame[column][0])):
        letters = []
        for item in range(0, len(pd_frame[column])):
            letters.append(pd_frame[column][item][letter])
        if all(x == letters[0] for x in letters) == True:
            conserved_region.append('*')
        else:
            conserved_region.append(' ')
    return conserved_region

def create_alignment_file_segment(pd_frame:'pd.DataFrame', column_name:'str', file_name:'str', segment:'str') -> Bio.Align.MultipleSeqAlignment:
    '''This function creates a fasta file suitable for alignment'''
    
    alignment_list = []
    min_legnth = pd_frame[column_name].map(lambda x: len(x[segment]) if len(x[segment]) > 0 else 10**5).min()
    
    for item in pd_frame['Record id']:
        if len(pd_frame.loc[(pd_frame['Record id'] == item)][column_name].item()[segment]) > 1:
            species_location = item.find('.')
            alignment_list.append(SeqRecord(Seq(pd_frame.loc[(pd_frame['Record id'] == item)][column_name].item()[segment][1:min_legnth], generic_dna), 
            id = item,
            description = str(item[0] + '.' + item[species_location + 1:])))
        else:
            pass
        
    writen_alignment = MultipleSeqAlignment(alignment_list)
    AlignIO.write(writen_alignment, str(file_name +'.faa'), "fasta")
    return writen_alignment  
