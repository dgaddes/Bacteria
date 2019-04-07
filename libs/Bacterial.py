#Bacterial
import os
import glob
import numpy as np
import pandas as pd
import Bio
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC

def get_species():
    os.chdir('/home/dave/Documents/Bacteria/16S') #$echo $PWD to determine the path name
    bacterial_species = []
    all_seq_records = []

    for file in glob.glob("*.fna"):
        bacterial_species.append(file)

    for file in bacterial_species:
        for seq_record in SeqIO.parse(file, "fasta"):
            all_seq_records.append([str(seq_record.description.split()[1]) + "." + str(seq_record.description.split()[2]),
            seq_record.id,
            str(seq_record.seq),
            len(seq_record)])
    return all_seq_records

def get_species_23S():
    os.chdir('/home/dave/Documents/Bacteria/23S') #$echo $PWD to determine the path name
    bacterial_species = []
    all_seq_records = []

    for file in glob.glob("*.fna"):
        bacterial_species.append(file)

    for file in bacterial_species:
        for seq_record in SeqIO.parse(file, "fasta"):
            all_seq_records.append([str(seq_record.description.split()[1]) + "." + str(seq_record.description.split()[2]),
            seq_record.id,
            str(seq_record.seq),
            len(seq_record)])
    return all_seq_records

def compare(amplicon_a, amplicon_b, shift):
    if len(amplicon_a) >= len(amplicon_b):
        reference_a = amplicon_a
        reference_b = amplicon_b
    else:
        reference_b = amplicon_a
        reference_a = amplicon_b
    comparision = str()
    for i in range(0, len(reference_b)):
        if reference_a[i+shift] == reference_b[i]:
            comparision = comparision + str('*')
        else:
            comparision = comparision + str('[' + reference_a[i] + '/' + reference_b[i] + ']')
    for i in range(0, len(reference_a) - len(reference_b)):
        comparision = comparision + str('[' + reference_a[i] + '/' + '$' + ']')
    print(comparision)

def compare_letters(amplicon_a, amplicon_b, shift):
    if len(amplicon_a) >= len(amplicon_b):
        reference_a = amplicon_a
        reference_b = amplicon_b
    else:
        reference_b = amplicon_a
        reference_a = amplicon_b
    comparision = str()
    for i in range(0, len(reference_b)):
        if reference_a[i+shift] == reference_b[i]:
            comparision = comparision + str(reference_a[i])
        else:
            comparision = comparision + str('[' + reference_a[i] + '/' + reference_b[i] + ']')
    for i in range(0, len(reference_a) - len(reference_b)):
        comparision = comparision + str('[' + reference_a[i] + '/' + '$' + ']')
    print(comparision)

def compare_fragments(amplicon_a, amplicon_b, basepair):
    fragment_a = []
    fragment_b = []
    if len(amplicon_a) >= len(amplicon_b):
        reference_a = amplicon_a
        reference_b = amplicon_b
    else:
        reference_b = amplicon_a
        reference_a = amplicon_b
    for i in range(0, len(reference_a)):
        if reference_a[i] == basepair:
            fragment_a.append(1)
        else:
            fragment_a.append(0)
    for i in range(0, len(reference_b)):
        if reference_b[i] == basepair:
            fragment_b.append(1)
        else:
            fragment_b.append(0)
    for i in range(len(reference_b), len(reference_a)):
    	fragment_b.append(0)


    return	fragment_a, fragment_b