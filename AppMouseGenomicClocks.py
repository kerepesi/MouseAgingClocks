import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr, ttest_rel, ttest_1samp
from scipy import stats
import sys

# Notes:
# - Example for usage: python AppMouseGenomicClocks.py FT.csv Meta.csv
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
Meta=pd.read_csv(IN2_name,index_col=[0])

## Multi-tissue epigenetic aging clock (Meer 2018) 
WLMTClock=pd.read_excel('ClockData/elife-40675-supp3-v2.xlsx', sheet_name='Whole lifespan multi-tissue', nrows=435)
WLMTClock['Pos'] = [ str(WLMTClock['Chromosome'][ind])+'_'+str(WLMTClock['Position'][ind]) 
                    for ind in WLMTClock.index ]
WLMTClock=WLMTClock.set_index('Pos')
clock_mean_mean=Clock_mean_mean(FT,WLMTClock.index)
Pred_dict={}
Overlap_dict={}
for Name in FT.columns:
    FT_app_dn=(FT*100)[[Name]].dropna()
    FT_app_reind_clock=FT_app_dn.reindex(WLMTClock.index)
    FT_app_reind_clock_fillna=FT_app_reind_clock.fillna(clock_mean_mean).T
    Pred=((FT_app_reind_clock_fillna)*WLMTClock['Weight']).sum(axis=1)+234.64
    Pred=Pred[0]/30.5
    n_overlapped_sites=len(WLMTClock.index.intersection(FT_app_dn.index))
    Pred_dict[Name]=Pred
    Overlap_dict[Name]=n_overlapped_sites
clock='Meer multi-tissue'
value_name='DNAm age ('+clock+')'
Meta[value_name]=pd.Series(Pred_dict)
Meta['Overlap '+clock]=pd.Series(Overlap_dict)

## Petkovich blood clock
BloodClock=pd.read_excel('ClockData/elife-40675-supp3-v2.xlsx', sheet_name='Blood', nrows=90)
d={}
for ind in BloodClock.index:
    pos=str(BloodClock['Chromosome'][ind])+'_'+str(BloodClock['Position'][ind]) 
    d[ind]=pos
BloodClock['Pos']=pd.Series(d)
BloodClock=BloodClock.set_index('Pos')
clock_mean_mean=Clock_mean_mean(FT,BloodClock.index)
Pred_dict={}
Overlap_dict={}
for Name in FT.columns:
    FT_app_dn=FT[[Name]].dropna()
    FT_app_reind_clock=FT_app_dn.reindex(BloodClock['Weight'].index)
    FT_app_reind_clock_fillna=FT_app_reind_clock.fillna(clock_mean_mean).T
    MSc=((FT_app_reind_clock_fillna)*BloodClock['Weight']).sum(axis=1)
    MSc=MSc[0]
    a = 0.1666
    b = 0.4185 
    c = - 1.712
    Age = ((MSc - c)/a)**(1/b)
    Pred=Age/30.5
    n_overlapped_sites=len(BloodClock['Weight'].index.intersection(FT_app_dn.index))
    Pred_dict[Name]=Pred
    Overlap_dict[Name]=n_overlapped_sites
clock='Petkovich blood'
value_name='DNAm age ('+clock+')'
Meta[value_name]=pd.Series(Pred_dict)
Meta['Overlap '+clock]=pd.Series(Overlap_dict)

## Stubbs multi-tissue clock
YOMT=pd.read_excel('ClockData/elife-40675-supp3-v2.xlsx', sheet_name='Young age multi-tissue', nrows=329)
d={}
for ind in YOMT.index:
    pos=str(YOMT['Chromosome'][ind])+'_'+str(YOMT['Position'][ind]) 
    d[ind]=pos
YOMT['Pos']=pd.Series(d)
YOMT=YOMT.set_index('Pos')
Train=pd.read_csv('ClockData/TrainingData_Babraham_Reizel_Cannon.txt', sep='\t',index_col=[0])
Train['Pos']=['chr'+ind.replace(':','_') for ind in Train.index]
Train=Train.set_index(['Pos'])
Train_reind=Train.reindex(YOMT.index).T
clock_mean_mean=Clock_mean_mean(FT,YOMT.index)
Pred_dict={}
Overlap_dict={}
for Name in FT.columns:
    FT_app_dn=FT[[Name]].dropna()
    FT_app_reind_clock=FT_app_dn.reindex(YOMT.index)
    FT_app_reind_clock_fillna=FT_app_reind_clock.fillna(clock_mean_mean).T
    FT_app_reind_clock_fillna_norm=(FT_app_reind_clock_fillna - Train_reind.median()) / Train_reind.std()
    MSc=((FT_app_reind_clock_fillna_norm)*YOMT['Weight']).sum(axis=1)
    MSc=MSc[0]
    Age = np.exp (0.1207 * (MSc**2) + 1.2424 * MSc + 2.5440) - 3
    Pred=Age*(7/30.5)
    n_overlapped_sites=len(YOMT.index.intersection(FT_app_dn.index))
    Pred_dict[Name]=Pred
    Overlap_dict[Name]=n_overlapped_sites
clock='Stubbs multi-tissue'
value_name='DNAm age ('+clock+')'
Meta[value_name]=pd.Series(Pred_dict)
Meta['Overlap '+clock]=pd.Series(Overlap_dict)

## Thompson multi-tissue EN
Weight={}
f=open('ClockData/Thompson2018-ElasticNet_aging_clock.txt')
i=-1
for l in f:
    i+=1
    sl=l.strip().split('\t')
    if i==2:
        Intercept=float(sl[2])
    elif i>2:
        Pos=sl[0]+'_'+sl[1]
        Weight[Pos]=float(sl[2])
Thompson_EN=pd.Series(Weight)
clock_mean_mean=Clock_mean_mean(FT,Thompson_EN.index)
Pred_dict={}
Overlap_dict={}
for Name in FT.columns:
    FT_app_dn=FT[[Name]].dropna()
    FT_app_reind_clock=FT_app_dn.reindex(Thompson_EN.index)
    FT_app_reind_clock_fillna=FT_app_reind_clock.fillna(clock_mean_mean).T
    Pred=(FT_app_reind_clock_fillna.values*Thompson_EN.values).sum(axis=1)+30.3172
    Pred=Pred[0]   
    n_overlapped_sites=len(Thompson_EN.index.intersection(FT_app_dn.index))
    Pred_dict[Name]=Pred
    Overlap_dict[Name]=n_overlapped_sites
clock='Thompson multi-tissue EN'
value_name='DNAm age ('+clock+')'
Meta[value_name]=pd.Series(Pred_dict)
Meta['Overlap '+clock]=pd.Series(Overlap_dict)

## Out
OUT_name=IN1_name+'-'+sys.argv[0].split('/')[-1]+'.csv'
print(OUT_name)
Meta.to_csv(OUT_name)
