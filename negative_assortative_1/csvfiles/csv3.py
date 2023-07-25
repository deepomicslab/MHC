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
from xml.etree import ElementTree as ET
from statsmodels.stats.multitest import multipletests





all_mean_arr=-0.017462454026226405
all_std_arr=0.006996042091888453
mean_dict_different_race={'ACB': -0.021129940274677197,
 'ASW': -0.020541364517036285,
 'BEB': -0.01690678177281044,
 'CDX': -0.01520099465128729,
 'CEU': -0.016533114461341817,
 'CHB': -0.015328450587149127,
 'CHS': -0.01553630769670323,
 'CLM': -0.017175762377800757,
 'ESN': -0.021289481956977087,
 'FIN': -0.016220012008386193,
 'GBR': -0.01616972796402513,
 'GIH': -0.016647828395962485,
 'GWD': -0.02127146261125899,
 'IBS': -0.016612028858784364,
 'ITU': -0.016483044568386403,
 'JPT': -0.015297644767223825,
 'KHV': -0.015510621661266998,
 'LWK': -0.020977118301868332,
 'MSL': -0.021221971390719675,
 'MXL': -0.016358508148574633,
 'PEL': -0.015192057641423692,
 'PJL': -0.016650973867877568,
 'PUR': -0.01772660210939756,
 'STU': -0.016488613724290634,
 'TSI': -0.016378152064038587,
 'YRI': -0.021395734138643974}
std_dict_different_race={'ACB': 0.007866911384050993,
 'ASW': 0.007845222560823078,
 'BEB': 0.007884710711424594,
 'CDX': 0.007874313523803581,
 'CEU': 0.007927585629080875,
 'CHB': 0.007946756589325447,
 'CHS': 0.00794164324551959,
 'CLM': 0.007812970518334978,
 'ESN': 0.00785082325216948,
 'FIN': 0.007906254940044924,
 'GBR': 0.007932404777081037,
 'GIH': 0.007873289737229847,
 'GWD': 0.007853536883558044,
 'IBS': 0.007919548085211019,
 'ITU': 0.007768497466178242,
 'JPT': 0.00796914647357699,
 'KHV': 0.007990414654274063,
 'LWK': 0.007797858211516251,
 'MSL': 0.007864701129783158,
 'MXL': 0.007798067190318037,
 'PEL': 0.007722619015518075,
 'PJL': 0.007712973460516593,
 'PUR': 0.007785891217808583,
 'STU': 0.0077724707332740475,
 'TSI': 0.00791282517960449,
 'YRI': 0.007915268495722334}

df_csv2=pd.read_csv("/data2/wangxuedong/mhc_test_data/similarity_region_pythonfiles/csvfiles/csv2_first_with_repeat_from_ipynb.csv").values.tolist()
#get indices of minimum value
def getmaxindices(arr,num):
    max_values_and_indices = sorted(enumerate(arr), key=lambda x: x[1])[:num]
    max_indices = [index for index, value in max_values_and_indices]    
    return max_indices

allrecords=[]
#in csv2
for record in df_csv2:
    if type(record[6])!=float:
        print("record",record)
        print("record6",record[6])
        if record[1]=="All": 
            if record[2]==0:#no mhc csvfiles     all+nomhc 
                print("this is all + nomhc.")                  
                position_of_gene=record[6]
                chr=position_of_gene.split(":")[0]
                start_position=int(position_of_gene.split(":")[1].split("-")[0])

                end_position=int(position_of_gene.split(":")[1].split("-")[1])
                filenames=os.listdir("/data2/wangxuedong/mhc_test_data/similar_region_csvfiles/")
                for filename in filenames:
                    if filename.startswith("pos_prob_split_"+chr+"_"):
                        df_split=pd.read_csv("/data2/wangxuedong/mhc_test_data/similar_region_csvfiles/"+filename)  

                        allpos=df_split[(df_split["Unnamed: 0"]>=start_position) & (df_split["Unnamed: 0"]<=end_position)]["Unnamed: 0"].to_list()
                        allrigion_P_top_1000=df_split[(df_split["Unnamed: 0"]>=start_position) & (df_split["Unnamed: 0"]<=end_position)]["P"].to_list()
                        if allrigion_P_top_1000:
                            indices_top_100=getmaxindices(allrigion_P_top_1000,100)
                            p_top_100=[allrigion_P_top_1000[i] for i in indices_top_100]
                            pos_top_100=[allpos[i] for i in indices_top_100]
                            for j in range(len(p_top_100)):
                                pvalue=round(stats.norm.cdf(p_top_100[j], all_mean_arr, all_std_arr),15)
                                onerecord=[]
                                onerecord+=record[1:8]
                                onerecord+=[pos_top_100[j],p_top_100[j],pvalue]
                                allrecords.append(onerecord)
                                

                        else:
                            continue

            else:#all+mhc
                position_of_gene=record[6]
                start_position=int(position_of_gene.split(":")[1].split("-")[0])
                end_position=int(position_of_gene.split(":")[1].split("-")[1])
                df_mhc=pd.read_csv("/data2/wangxuedong/mhc_test_data/similar_region_csvfiles/posmhc.csv")
                allpos=df_mhc[(df_mhc["Unnamed: 0"]>=start_position) & (df_mhc["Unnamed: 0"]<=end_position)]["Unnamed: 0"].to_list()
                allrigion_P_top_1000=df_mhc[(df_mhc["Unnamed: 0"]>=start_position) & (df_mhc["Unnamed: 0"]<=end_position)]["P"].to_list()
                if allrigion_P_top_1000:
                    indices_top_100=getmaxindices(allrigion_P_top_1000,100)
                    p_top_100=[allrigion_P_top_1000[i] for i in indices_top_100]
                    pos_top_100=[allpos[i] for i in indices_top_100]
                    for j in range(len(p_top_100)):
                        pvalue=round(stats.norm.cdf(p_top_100[j], all_mean_arr, all_std_arr),15)
                        onerecord=[]
                        onerecord+=record[1:8]
                        onerecord+=[pos_top_100[j],p_top_100[j],pvalue]
                        allrecords.append(onerecord)
                else:
                    continue

        else:#race!=all
            if record[2]==0:#no mhc race+nomhc
                position_of_gene=record[6]
                chr=position_of_gene.split(":")[0]
                start_position=int(position_of_gene.split(":")[1].split("-")[0])

                end_position=int(position_of_gene.split(":")[1].split("-")[1])
                filenames=os.listdir("/data2/wangxuedong/mhc_test_data/similar_region_race/")
                for filename in filenames:
                    if filename.startswith("split_"+chr+"_"):
                        df_split=pd.read_csv("/data2/wangxuedong/mhc_test_data/similar_region_race/"+filename)  
                        allpos=df_split[(df_split["Unnamed: 0"]>=start_position) & (df_split["Unnamed: 0"]<=end_position)]["Unnamed: 0"].to_list()
                        allrigion_P_top_1000=df_split[(df_split["Unnamed: 0"]>=start_position) & (df_split["Unnamed: 0"]<=end_position)][record[1]].to_list()
                        if allrigion_P_top_1000:
                            indices_top_100=getmaxindices(allrigion_P_top_1000,100)
                            p_top_100=[allrigion_P_top_1000[i] for i in indices_top_100]
                            pos_top_100=[allpos[i] for i in indices_top_100]
                            for j in range(len(p_top_100)):
                                pvalue=round(stats.norm.cdf(p_top_100[j],mean_dict_different_race[record[1]],std_dict_different_race[record[1]]),15)
                                onerecord=[]
                                onerecord+=record[1:8]
                                onerecord+=[pos_top_100[j],p_top_100[j],pvalue]
                                allrecords.append(onerecord)
                        else:
                            continue
                    else:
                        continue
            else:#ismhc different race
                position_of_gene=record[6]
                start_position=int(position_of_gene.split(":")[1].split("-")[0])
                end_position=int(position_of_gene.split(":")[1].split("-")[1])
                df_mhc=pd.read_csv("/data2/wangxuedong/mhc_test_data/similar_region_race/mhc.csv")
                allpos=df_mhc[(df_mhc["Unnamed: 0"]>=start_position) & (df_mhc["Unnamed: 0"]<=end_position)]["Unnamed: 0"].to_list()
                allrigion_P_top_1000=df_mhc[(df_mhc["Unnamed: 0"]>=start_position) & (df_mhc["Unnamed: 0"]<=end_position)][record[1]].to_list()
                if allrigion_P_top_1000:   
                    indices_top_100=getmaxindices(allrigion_P_top_1000,100)
                    p_top_100=[allrigion_P_top_1000[i] for i in indices_top_100]    
                    pos_top_100=[allpos[i] for i in indices_top_100]     
                    for j in range(len(p_top_100)):
                        pvalue=round(stats.norm.cdf(p_top_100[j],mean_dict_different_race[record[1]],std_dict_different_race[record[1]]),15)
                        onerecord=[]
                        onerecord+=record[1:8]
                        onerecord+=[pos_top_100[j],p_top_100[j],pvalue]
                        allrecords.append(onerecord)
                else:
                    continue
    else:
        print("record",record)
        print("record6",record[6])
        onerecord=[]
        onerecord+=record[1:8]
        onerecord+=["null","null","null"]
        allrecords.append(onerecord)


#generatedf
def generatedf(columns,allrecords):
    dictforDF=dict()
    for i in range(len(columns)):
        midarr=[]
        for record in allrecords:           
            midarr.append(record[i])
        dictforDF[columns[i]]=midarr
    dataframe(dictforDF).to_csv("/data2/wangxuedong/mhc_test_data/similarity_region_pythonfiles/csvfiles/csv3_top100_without_corr_frompy.csv")
    return dataframe(dictforDF)


dfnew=generatedf(columns=["Race","isMHC","Position of disassortative mating region","Gene ID","Gene name","Position of gene","is_complement","Position of SNP","P","p-value"],allrecords=allrecords)


pvalue_list=dfnew["p-value"].tolist()
corr_pvals = multipletests(pvalue_list, method='fdr_bh')[1]
dfnew = dfnew.assign(corr_pvals=corr_pvals)

dfnew.to_csv("/data2/wangxuedong/mhc_test_data/similarity_region_pythonfiles/csvfiles/csv3_top100_with_corr_frompy.csv")
              