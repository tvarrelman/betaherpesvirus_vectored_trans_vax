#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 17:39:08 2021

@author: tannervarrelman
"""

import pandas as pd 
import numpy as np
from statistics import mode

sims = pd.read_csv('/Users/tannervarrelman/Documents/MCMV_project_5_20/DurationPathData/Figure3Sim_6_28_21.csv')
# sims = pd.read_csv('/Users/tannervarrelman/Documents/MCMV_project_5_20/DurationPathData/Figure3Sim_2_24_21.csv')
# sims = pd.read_csv('/Users/tannervarrelman/Documents/MCMV_project_5_20/DurationPathData/Figure3Sim_reduced_eff_2_24_21.csv')
# sims = pd.read_csv('/Users/tannervarrelman/Documents/MCMV_project_5_20/DurationPathData/Figure3Sim_reduc_3_8_21.csv')
LASV_sub = sims[sims['Path.label']=='Lassa virus'].reset_index(drop=True)
LCMV_sub = sims[sims['Path.label']=='Lymphocytic choriomeningitis virus'].reset_index(drop=True)

new_LASV = pd.DataFrame(columns = ['minP', 'maxP', 'meanP', 'time',
                                   'ReducMin', 'ReducMax', 'ReducMean'])
new_LCMV = pd.DataFrame(columns = ['minP', 'maxP', 'meanP', 'time',
                                   'ReducMin', 'ReducMax', 'ReducMean'])

d = 0.002739726
b = 500*d
LASV_delt = 0.0457
LASV_beta = 0.072659585

LCMV_delt = 0.002739726
LCMV_beta = 0.006052883

SSP = b*((1/(d+LASV_delt))-(1/LASV_beta))*0.05
LCMV_SSP = b*((1/(d+LCMV_delt))-(1/LCMV_beta))*0.05

# time = LASV_sub[LASV_sub['P'] <= SSP]

reduc_list = []
for sim in LASV_sub['SimNum'].unique():
    sim_sub = LASV_sub[LASV_sub['SimNum']==sim].reset_index(drop=True)
    reduc_time = sim_sub[sim_sub['P'] <= SSP]['time'].tolist()
    # print(reduc_time)
    reduc_list.append(reduc_time[0])
min_LASV_time = np.min(reduc_list)
max_LASV_time = np.max(reduc_list)
mean_LASV_time = np.mean(reduc_list)

LCMV_reduc_list = []
for sim2 in LCMV_sub['SimNum'].unique():
    LCMV_sim_sub = LCMV_sub[LCMV_sub['SimNum']==sim2].reset_index(drop=True)
    LCMV_reduc_time = LCMV_sim_sub[LCMV_sim_sub['P'] <= LCMV_SSP]['time'].tolist()
    # print(reduc_time)
    LCMV_reduc_list.append(LCMV_reduc_time[0])
min_LCMV_time = np.min(LCMV_reduc_list)
max_LCMV_time = np.max(LCMV_reduc_list)
mean_LCMV_time = np.mean(LCMV_reduc_list)


for time in range(0,3000):
    LASV_t_sub = LASV_sub[LASV_sub['time']==time].reset_index(drop=True)
    LASVminP = np.min(LASV_t_sub['P'])
    LASVmaxP = np.max(LASV_t_sub['P'])
    LASVmeanP = np.mean(LASV_t_sub['P'])
    new_LASV = new_LASV.append({'minP':LASVminP,
                                'maxP':LASVmaxP,
                                'meanP':LASVmeanP,
                                'time':time,
                                'ReducMin':min_LASV_time,
                                'ReducMax':max_LASV_time,
                                'ReducMean':mean_LASV_time}, ignore_index=True)
    
    LCMV_t_sub = LCMV_sub[LCMV_sub['time']==time].reset_index(drop=True)
    LCMVminP = np.min(LCMV_t_sub['P'])
    LCMVmaxP = np.max(LCMV_t_sub['P'])
    LCMVmeanP = np.mean(LCMV_t_sub['P'])
    new_LCMV = new_LCMV.append({'minP':LCMVminP,
                                'maxP':LCMVmaxP,
                                'meanP':LCMVmeanP,
                                'time':time,
                                'ReducMin':min_LCMV_time,
                                'ReducMax':max_LCMV_time,
                                'ReducMean':mean_LCMV_time}, ignore_index=True)

new_LCMV.to_csv('/Users/tannervarrelman/Documents/MCMV_project_5_20/DurationPathData/LCMV_figure3_6_28_21.csv', header=True, index=False)
new_LASV.to_csv('/Users/tannervarrelman/Documents/MCMV_project_5_20/DurationPathData/LASV_figure3_6_28_21.csv', header=True, index=False)








