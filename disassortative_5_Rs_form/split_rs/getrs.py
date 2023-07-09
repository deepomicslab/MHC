#!/home/wangxuedong/app/miniconda3/envs/mhc/bin/python
import pysam
from pysam import VariantFile as vcf
import operator
from math import log2
import pandas as pd
from pandas import DataFrame as dataframe
import matplotlib.pyplot as plt
import numpy as np

from scipy.spatial.distance import pdist, squareform
import scipy
import  os
import os.path
import matplotlib.colors as mcolors
from scipy import stats
import csv
import requests
import json
from xml.etree import ElementTree as ET
from statsmodels.stats.multitest import multipletests
import re
import pickle
split_race_region={'ACB': {'split_chr2/xau': [[89762579, 89852967]],
  'split_chr5/xaa': [[676694, 846549]],
  'split_chr7/xan': [[56467013, 56713542]],
  'split_chr11/xaw': [[103206633, 103451621]],
  'split_chr14/xaa': [[19806333, 19976578]],
  'split_chr21/xaa': [[10605321, 10733724]]},
  'All':{'split_chr14/xaa': [[19806333, 19976578]],
 'split_chr7/xbc': [[124810550, 125070920]],
 'split_chr1/xak': [[45762038, 46100296]],
 'split_chr1/xbk': [[188243362, 188513984]]},
 'ASW': {'split_chr2/xau': [[89762579, 89852967]],
  'split_chr5/xaa': [[676694, 846549]],
  'split_chr9/xak': [[42796452, 42969270]],
  'split_chr10/xaj': [[38491573, 38687800], [38725237, 38832458]],
  'split_chr11/xaf': [[23272931, 23495433]],
  'split_chr13/xaa': [[18347994, 18534139]],
  'split_chr14/xak': [[66346187, 66752333]],
  'split_chr21/xaa': [[10605321, 10733724]]},
 'BEB': {'split_chr1/xbx': [[245939605, 246162003]],
  'split_chr2/xan': [[56405225, 56644106]],
  'split_chr3/xav': [[98046342, 98273693]],
  'split_chr3/xbf': [[145365797, 145612467]],
  'split_chr4/xaz': [[115247002, 115485969]],
  'split_chr13/xai': [[55102733, 55443491]],
  'split_chr14/xaa': [[19806333, 19976578]]},
 'CDX': {'split_chr1/xay': [[114816228, 115116610]],
  'split_chr2/xbd': [[130125757, 130339027]],
  'split_chr4/xba': [[119293026, 119539168]],
  'split_chr6/xao': [[64160994, 64523871]],
  'split_chr7/xaf': [[19734421, 19953255]],
  'split_chr9/xau': [[105652970, 105876427]],
  'split_chr10/xan': [[55664607, 55894887]],
  'split_chr12/xbc': [[131390693, 131619328]],
  'split_chr13/xai': [[55102733, 55354096]],
  'split_chr17/xag': [[31068982, 31384157]],
  'split_chr18/xan': [[65996302, 66233011]]},
 'CEU': {'split_chr1/xbe': [[158399387, 158641706]],
  'split_chr6/xac': [[8413699, 8649798]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr11/xaf': [[24934206, 25144712]],
  'split_chr12/xac': [[10962097, 11205540]],
  'split_chr17/xaf': [[21588491, 21823237]],
  'split_chr17/xaj': [[45620566, 46418024]]},
 'CHB': {'split_chr2/xbd': [[130125757, 130339027]],
  'split_chr4/xba': [[119293026, 119587470]],
  'split_chr4/xaz': [[115247002, 115485969]],
  'split_chr6/xao': [[64202650, 64523871]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr7/xaf': [[19775051, 19953255]],
  'split_chr10/xar': [[73060579, 73337762]],
  'split_chr11/xav': [[98072700, 98330585]],
  'split_chr22/xaa': [[16175432, 16345129]]},
 'CHS': {'split_chr3/xav': [[98046342, 98245327]],
  'split_chr4/xba': [[119293026, 119539168]],
  'split_chr5/xaa': [[676694, 846549]],
  'split_chr6/xao': [[64202650, 64523871]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr8/xaz': [[109058766, 109345554]],
  'split_chr14/xat': [[105509278, 105695533]]},
 'CLM': {'split_chr1/xbd': [[152162701, 152471006]],
  'split_chr3/xav': [[95059735, 95322404], [98046342, 98273693]],
  'split_chr3/xbj': [[163862679, 164224214]],
  'split_chr4/xba': [[119293026, 119587470]],
  'split_chr6/xab': [[7345120, 7558507]],
  'split_chr6/xar': [[76169005, 76460170]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr12/xan': [[61013177, 61213757]],
  'split_chr14/xaa': [[19806333, 19976578]],
  'split_chr14/xat': [[106520496, 106674883]]},
 'ESN': {'split_chr2/xau': [[89581761, 89852967]],
  'split_chr5/xaa': [[676694, 846549]],
  'split_chr5/xai': [[34311870, 34562837]],
  'split_chr6/xbf': [[145443018, 145762353]],
  'split_chr7/xap': [[67044394, 67267582]],
  'split_chr14/xan': [[80889782, 81152442]],
  'split_chr21/xaa': [[10605321, 10733724]]},
 'FIN': {'split_chr3/xbj': [[163862679, 164147845]],
  'split_chr4/xaw': [[97892649, 98142264]],
  'split_chr4/xbh': [[151414350, 151728888]],
  'split_chr6/xat': [[85153819, 85387984]],
  'split_chr6/xbf': [[145485121, 145762353]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr8/xaz': [[109010721, 109345554]],
  'split_chr9/xav': [[112161299, 112408537]],
  'split_chr11/xbb': [[124171937, 124450010]],
  'split_chr18/xac': [[12017790, 12244103]]},
 'GBR': {'split_chr1/xak': [[45628366, 46147211]],
  'split_chr2/xan': [[56405225, 56644106]],
  'split_chr2/xas': [[78325225, 78577282]],
  'split_chr3/xbf': [[145365797, 145612467]],
  'split_chr3/xbi': [[158122703, 158435520]],
  'split_chr4/xba': [[119243365, 119539168]],
  'split_chr6/xac': [[8413699, 8649798]],
  'split_chr7/xbc': [[124760998, 125070920]],
  'split_chr11/xbb': [[124171937, 124450010]],
  'split_chr12/xac': [[10962097, 11205540]],
  'split_chr17/xaf': [[21654230, 21823237]],
  'split_chr17/xaj': [[45587933, 46418024]]},
 'GIH': {'split_chr3/xav': [[98046342, 98276757]],
  'split_chr4/xbd': [[132580391, 132816794]],
  'split_chr6/xac': [[8413699, 8649798]],
  'split_chr11/xah': [[34650776, 34905225]],
  'split_chr14/xaa': [[19806333, 20015421]],
  'split_chr18/xac': [[12017790, 12244103]]},
 'GWD': {'split_chr1/xak': [[45803471, 46100296]],
  'split_chr2/xai': [[38264369, 38464894]],
  'split_chr2/xau': [[89762579, 89852967]],
  'split_chr5/xaa': [[676694, 846549]],
  'split_chr13/xaa': [[18347994, 18534139]],
  'split_chr21/xaa': [[10511846, 10763989]]},
 'IBS': {'split_chr3/xbf': [[145365797, 145612467]],
  'split_chr6/xbf': [[145485121, 145762353]],
  'split_chr12/xac': [[10962097, 11205540]],
  'split_chr12/xas': [[83579779, 83843102]],
  'split_chr17/xaf': [[21619971, 21823237]],
  'split_chr17/xaj': [[45836662, 46418024]],
  'split_chr18/xac': [[12017790, 12244103]],
  'split_chr22/xaa': [[16131708, 16345129]]},
 'ITU': {'split_chr5/xas': [[84862293, 85266031]],
  'split_chr10/xan': [[55664607, 55894887]],
  'split_chr12/xac': [[10962097, 11162840]],
  'split_chr14/xaa': [[19806333, 20047342]],
  'split_chr17/xak': [[52575391, 52889013]]},
 'JPT': {'split_chr3/xav': [[98046342, 98245327]],
  'split_chr5/xbc': [[131688028, 132052794]],
  'split_chr5/xaa': [[676694, 825167]],
  'split_chr13/xaa': [[18347994, 18534139]],
  'split_chr13/xap': [[85638733, 85884384]],
  'split_chr14/xad': [[31681845, 31909184]],
  'split_chr22/xaa': [[16175432, 16345129]]},
 'KHV': {'split_chr2/xbq': [[192771748, 193054217]],
  'split_chr3/xav': [[98046342, 98245327]],
  'split_chr3/xbj': [[163946649, 164185047]],
  'split_chr4/xba': [[119293026, 119539168]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr14/xak': [[66298204, 66793088]]},
 'LWK': {'split_chr2/xai': [[38235897, 38464894]],
  'split_chr2/xau': [[89762579, 89852967]],
  'split_chr5/xaa': [[676694, 846549]],
  'split_chr13/xaa': [[18347994, 18534139]],
  'split_chr14/xak': [[66346187, 66752333]],
  'split_chr21/xaa': [[10605321, 10733724]]},
 'MSL': {'split_chr1/xak': [[45803471, 46100296]],
  'split_chr2/xau': [[89762579, 89852967]],
  'split_chr5/xay': [[113365150, 113639918]],
  'split_chr7/xaf': [[19775051, 19953255]],
  'split_chr7/xap': [[67044394, 67267582]],
  'split_chr9/xak': [[42796452, 42969270]],
  'split_chr11/xaf': [[23272931, 23495433]],
  'split_chr11/xaw': [[103206633, 103451621]],
  'split_chr13/xaa': [[18347994, 18534139]],
  'split_chr14/xat': [[105569851, 105695533]],
  'split_chr21/xaa': [[10605321, 10733724]]},
 'MXL': {'split_chr3/xav': [[98046342, 98217489]],
  'split_chr3/xaq': [[74901608, 75192541]],
  'split_chr5/xad': [[11999758, 12263260]],
  'split_chr6/xao': [[64160994, 64404669]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr11/xaf': [[23200593, 23461651]],
  'split_chr12/xad': [[14846240, 15093187]],
  'split_chr12/xac': [[10962097, 11162840]],
  'split_chr12/xap': [[73338587, 73586502]],
  'split_chr12/xas': [[83579779, 83843102]],
  'split_chr17/xaf': [[21654230, 21823237]]},
 'PEL': {'split_chr1/xak': [[45628366, 46147211]],
  'split_chr2/xbo': [[186062200, 186331668]],
  'split_chr3/xav': [[98046342, 98217489]],
  'split_chr3/xba': [[120872663, 121116904]],
  'split_chr4/xaq': [[69225334, 69589122]],
  'split_chr6/xab': [[7345120, 7695065]],
  'split_chr8/xba': [[112137591, 112486758]],
  'split_chr10/xaq': [[68153912, 68481883]],
  'split_chr11/xaf': [[23168664, 23461651]],
  'split_chr11/xam': [[58630320, 58868016]],
  'split_chr12/xac': [[10962097, 11162840]],
  'split_chr12/xas': [[83579779, 83843102]],
  'split_chr14/xaa': [[19806333, 19976578]],
  'split_chr16/xab': [[5659954, 5828695]]},
 'PJL': {'split_chr12/xac': [[10962097, 11162840]],
  'split_chr14/xac': [[28519266, 28761965]],
  'split_chr14/xaa': [[19806333, 20078222]],
  'split_chr17/xak': [[52575391, 52861780]]},
 'PUR': {'split_chr3/xav': [[98003499, 98273693]],
  'split_chr3/xbj': [[163819938, 164063747]],
  'split_chr12/xac': [[10962097, 11205540]],
  'split_chr12/xan': [[61013177, 61213757]],
  'split_chr12/xas': [[83689105, 83928389]],
  'split_chr17/xaf': [[21654230, 21823237]]},
 'STU': {'split_chr3/xav': [[98046342, 98273693]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr13/xaa': [[18347994, 18534139]],
  'split_chr13/xai': [[55102733, 55443491]],
  'split_chr14/xaa': [[19806333, 19976578]],
  'split_chr16/xaj': [[48512974, 48799209]],
  'split_chr19/xam': [[51801039, 52098595]]},
 'TSI': {'split_chr4/xaw': [[97892649, 98189302]],
  'split_chr4/xba': [[119293026, 119539168]],
  'split_chr7/xag': [[23951217, 24204307]],
  'split_chr7/xbc': [[124810550, 125070920]],
  'split_chr9/xai': [[31957257, 32172669]],
  'split_chr12/xac': [[10962097, 11162840]],
  'split_chr14/xaa': [[19806333, 19976578]],
  'split_chr17/xaf': [[21654230, 21823237]],
  'split_chr17/xaj': [[45836662, 46544765]]},
 'YRI': {'split_chr2/xau': [[89762579, 89852967]],
  'split_chr5/xaa': [[676694, 846549]],
  'split_chr7/xan': [[56467013, 56713542]],
  'split_chr7/xap': [[67044394, 67267582]],
  'split_chr9/xak': [[42796452, 42969270]],
  'split_chr13/xaa': [[18347994, 18534139]],
  'split_chr16/xao': [[72288822, 72684953]],
  'split_chr21/xaa': [[10605321, 10733724]]}}



result_dict=dict()
for race, region_ in split_race_region.items():
    for split,regions in region_.items():
        for region in regions:
            result_dict[tuple(region)]=split

df=pd.read_csv("/data2/wangxuedong/mhc_test_data/python/getRS/df_sorted_0513_top100_csv3.csv")
df_split= df[df['isMHC'] == 0]



def getwhichsplit(resultdict,position,chr):
    allsplit=[]
    for digitregion,split_ in resultdict.items():
        if position>=digitregion[0] and position<=digitregion[1]:
            s=split_.split("/")[0].split("_")[-1]
            print(s)
            if chr==s:
                allsplit.append(split_)
        else:
            continue
    return list(set(allsplit))


with open('/data2/wangxuedong/mhc_test_data/getRS_0616/com_split.pkl', 'rb') as f:
    allrecords = pickle.load(f)

def extract_numbers(input_string):
    return re.findall(r'\d+', input_string)

outputs=[]

for arr in allrecords:
    whichsplit=getwhichsplit(result_dict,arr[1],arr[0])
    if len(whichsplit)==1:
        filepath="/data2/wangxuedong/mhc_test_data/"+whichsplit[0].split("/")[0]+"/"+whichsplit[0].split("/")[1]+".vcf.gz"
        print(filepath)
        vcf_in=vcf(filepath)
        for rec in vcf_in.fetch():
            if rec.pos==arr[1]:
                onerecord=["chromosome"]
                onerecord+=[extract_numbers(arr[0])[0],rec.pos,rec.ref,rec.alts[0],1]
                outputs.append(onerecord)
    else:
        print("with2split",arr)

def generatedf(columns,allrecords):
    dictforDF=dict()
    for i in range(len(columns)):
        midarr=[]
        for record in allrecords:           
            midarr.append(record[i])
        dictforDF[columns[i]]=midarr
    dataframe(dictforDF).to_csv("/data2/wangxuedong/mhc_test_data/getRS_0616/forms/COM_SPLIT_PART0.csv")
    return dataframe(dictforDF)


df2=generatedf(columns=["chrom","chr","pos","ref","alts[0]","not know"],allrecords=outputs)