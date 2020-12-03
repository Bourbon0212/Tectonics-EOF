# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 21:47:48 2020

@author: dx788
"""

import os
import math
import pandas as pd
import numpy as np

### 更換天線盤的測站

antenna = 'E:/GitHub/Tectonics-EOF/antenna.txt'
antenna = pd.read_csv(antenna, sep=r'\s*(?:\t|\s)\s*', header=None)
antenna = pd.DataFrame(np.concatenate((antenna.iloc[:,0:3].values, antenna.iloc[:,[0,3,4]].values), axis=0))
antenna[3] = antenna[1] + (antenna[2] - 0.5)/366

year_filter = (antenna[3] >= 2003.77459) & (antenna[3] <= 2004.10519)
antenna_f = antenna[year_filter].drop_duplicates()

### 讀檔

file_path = [x[:4] for x in os.listdir('data_raw')]
main_sta = list(np.setdiff1d(file_path, antenna_f[0]))

def least_square(start, end, event, component):
    beta = {}
    error = []
    for sta in main_sta:
        try:
            df = pd.read_csv('E:/GitHub/Tectonics-EOF/data_raw/' + sta + '.COR', header=None,
                            sep=r'\s*(?:\t|\s)\s*', engine='python',
                            names=['time','lat','lon','hgt','E','N','U','X'])
            # 研究取樣時間
            time_filter = (df.time >= start) & (df.time <= end)
            df_C = df.loc[time_filter,['time', component]].reset_index()
            
            # start_julian = math.ceil((start - math.floor(start)) * 366 + 0.5)
            # end_julian = math.ceil((end - math.floor(end)) * 366 + 0.5)
            
            # 常數
            constant = pd.DataFrame(np.repeat(1, len(df_C)))
            
            # 階梯函數
            step0_filter = (df.time >= start) & (df.time < event)
            step1_filter = (df.time >= event) & (df.time <= end)
            step = pd.concat([pd.DataFrame(np.repeat(0, sum(step0_filter))), pd.DataFrame(np.repeat(1, sum(step1_filter)))], axis=0, ignore_index=True)
            
            # 整理
            df_C = pd.concat([constant, step, df_C], axis=1, ignore_index=True).iloc[:,[0,1,3,4]]
            df_C.columns = ['cons','step','time','U']
            G = df_C.iloc[:,[0,1,2]].values; d = df_C.iloc[:,3].values
            
            # beta 係數
            m = np.linalg.inv(G.T@G)@(G.T@d) # (X'X)^(-1)X'y
            beta[sta] = m
        except:
            error.append(sta)
    return beta, error

start = 2003.77459; end = 2004.10519; event = 2003.93852
result = least_square(start, end, event, "U") 
        