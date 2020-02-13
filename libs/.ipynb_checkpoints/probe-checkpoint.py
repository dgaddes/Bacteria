#probe
'''
__author__ = "David Gaddes"
__Copyright__ "Copyright September 2019, David Gaddes"
__License__ = "GPL"
__email__ "dgaddes@protonmail.com"  
'''

import os, glob
import numpy as np
import pandas as pd
import Bio
from Bio.Seq import MutableSeq, Seq
from Bio import SeqIO
from typing import Tuple

offset = 100
variable_regions = ((0,0), (69-offset,99-offset), (137-offset,242-offset), (433-offset,497-offset), 
                    (576-offset,682-offset), (822-offset,879-offset), (986-offset,1043-offset),
                     (1117-offset,1173-offset),(1243-offset,1294-offset), (1435-offset,1465-offset))
conserved_regions_adjusted= [(0, 75), 
                             (112, 156),
                             (270, 471),
                             (536, 618), 
                             (725, 868), 
                             (928, 1038), 
                             (1096, 1176), 
                             (1233, 1305),
                             (1358, 1501)]
variable_regions_adjusted = [(75, 111), 
                             (156, 269), 
                             (471, 535), 
                             (618, 724), 
                             (868, 927), 
                             (1038, 1095), 
                             (1176, 1232), 
                             (1305, 1357), 
                             (1501, 1546)]

def splice_variable_regions(pd_frame:'pd.DataFrame', column:'str', start_reg:'str', end_reg:'str') -> pd.DataFrame:
    conserved_region = []
    variable_region = []
    
    for item in pd_frame['Record id']:
        sequence = pd_frame.loc[(pd_frame['Record id'] == item)][column].item()
        the_tuple = ()
        try:
            for region in range(start_reg, end_reg-start_reg+1):
                the_tuple = the_tuple + (sequence[variable_regions[region][0]:variable_regions[region][1]], )
            variable_region.append([pd_frame.loc[(pd_frame['Record id'] == item)]['Species'].item(),
                                   item,
                                   the_tuple])
        except:
            for region in range(start_reg, end_reg-start_reg+1):
                the_tuple = the_tuple + ('NaN', )

    columns = ['Species', 'Record id', 'Variable Regions',]
    returned_pd_frame = pd.DataFrame(variable_region, columns=columns)
    return returned_pd_frame