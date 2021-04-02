import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr, ttest_rel, ttest_1samp
from scipy import stats
import sys

# Notes:
# - Example for usage: python AppMouserDNAClocks.py FT.csv Meta.csv
# - Where FT.csv is a transposed feature table (rows: CpG sites, columns: samples)
#   and Meta.csv is a ',' separated metadata file (rows: samples, columns: attributes).
#   The list of samples should be the same as in FT.csv
# - Missing values are filled by the mean of the mean values of the clock sites of the samples

def Clock_mean_mean(FT, Clock_ind):
    l=[]
    for Name in FT.columns:
        FT_app_dn=FT[[Name]].dropna()
        FT_app_reind_clock=FT_app_dn.reindex(Clock_ind)
        l.append(FT_app_reind_clock.mean()[0])
    clock_mean_mean=np.mean(l)
    return clock_mean_mean

## Input
IN1_name=sys.argv[1]
IN2_name=sys.argv[2]
FT=pd.read_csv(IN1_name,index_col=[0])
Meta=pd.read_csv(IN2_name,index_col=[0],sep=',')

## Blood rDNA clock
WLClock=pd.read_excel('ClockData/Supplemental_Table_S2.xlsx',skiprows=2)
Intercept=WLClock['Weight'][0]
Weights=WLClock.drop(0)
l=[]
for ind in Weights.index:
    WL_pos=int(Weights['Unnamed: 0'][ind].split(',')[1])
    if WL_pos >= 1:
        original_pos=WL_pos
    else:
        original_pos=45306+WL_pos+1
    l.append(original_pos)
Weights['Original_pos']=l
Weights=Weights.set_index('Original_pos')
clock_mean_mean=Clock_mean_mean(FT,Weights.index)
Pred_dict={}
Overlap_dict={}
for Name in FT.columns:
    FT_app_dn=FT[[Name]].dropna()
    FT_app_reind_clock=FT_app_dn.reindex(Weights.index)
    FT_app_reind_clock_fillna=FT_app_reind_clock.fillna(clock_mean_mean).T
    Pred=(FT_app_reind_clock_fillna*Weights['Weight'].values).sum(axis=1)+Intercept
    Pred=np.exp(1)**Pred[0]
    n_overlapped_sites=len(Weights.index.intersection(FT_app_dn.index))
    Pred_dict[Name]=Pred
    Overlap_dict[Name]=n_overlapped_sites
clock='Blood rDNA'
value_name='rDNAm age ('+clock+')'
Meta[value_name]=pd.Series(Pred_dict)
Meta['Overlap '+clock]=pd.Series(Overlap_dict)

## Multi-tissue rDNA clock
MultiTissuerDNA=pd.read_csv('ClockData/MultiTissue_rDNA_clock-Overlap_ThompsonPetkovich-woMuscle-onlyC57BL-Lambda_opt-07112020.csv', header=None, index_col=[0])
Intercept = MultiTissuerDNA[1][0]
Weights=MultiTissuerDNA.drop('Intercept')
Weights.index=Weights.index.astype(int)
clock_mean_mean=Clock_mean_mean(FT,Weights.index)
Pred_dict={}
Overlap_dict={}
for Name in FT.columns:
    FT_app_dn=FT[[Name]].dropna()
    FT_app_reind_clock=FT_app_dn.reindex(Weights.index)
    FT_app_reind_clock_fillna=FT_app_reind_clock.fillna(clock_mean_mean).T
    Pred=(FT_app_reind_clock_fillna*(Weights[1].values)).sum(axis=1)+Intercept
    n_overlapped_sites=len(Weights.index.intersection(FT_app_dn.index))
    Pred_dict[Name]=Pred[0]
    Overlap_dict[Name]=n_overlapped_sites
clock='Multi-tissue rDNA'
value_name='rDNAm age ('+clock+')'
Meta[value_name]=pd.Series(Pred_dict)
Meta['Overlap '+clock]=pd.Series(Overlap_dict)

## Out
OUT_name=IN1_name+'-'+sys.argv[0].split('/')[-1]+'.csv'
print(OUT_name)
Meta.to_csv(OUT_name)
