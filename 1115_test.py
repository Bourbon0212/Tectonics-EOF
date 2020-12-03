# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 23:24:45 2020

@author: dx788
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 2003-12-10 (2003.93852 , 2003 344) 成功地震 c(121.4, 23.07)
# (2003.77459~2004.10519) 前後60天

df = pd.read_csv('E:/GitHub/Tectonics-EOF/data_raw/ANBU.COR', header=None,
                    sep=r'\s*(?:\t|\s)\s*', engine='python',
                    names=['time','lat','lon','hgt','E','N','U','X'])

time_filter = (df.time >= 2003.77459) & (df.time <= 2004.10519)
df_U = df.loc[time_filter,['time', 'U']].reset_index()

constant = pd.DataFrame(np.repeat(1, len(df_U)))
step = pd.concat([pd.DataFrame(np.repeat(0, 55)), pd.DataFrame(np.repeat(1, 59))], axis=0, ignore_index=True)

df_U = pd.concat([constant, step, df_U], axis=1, ignore_index=True).iloc[:,[0,1,3,4]]
df_U.columns = ['cons','step','time','U']

G = df_U.iloc[:,[0,1,2]].values; d = df_U.iloc[:,3].values
m = np.linalg.inv(G.T@G)@(G.T@d) # (X'X)^(-1)X'y

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')

# Generate a scatter plot
_ = df_U.plot(kind='scatter', x='time', y='U')
_ = abline(m[2],m[0])
_ = abline(m[2],m[0]+m[1])
plt.show()