# -*- coding: utf-8 -*-
"""
@author: slawler
"""

import pandas as pd
import requests
import json
from datetime import datetime
from collections import OrderedDict, defaultdict
import string
import requests
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
import os
from scipy.integrate import trapz, cumtrapz, simps
import numpy as np
from IPython.display import Markdown, display

def printbold(string):
    display(Markdown('**'+string+'**'))

def Get_USGS_Peaks(gage):  
    url = 'https://nwis.waterdata.usgs.gov/ny/nwis/peak?site_no={}&agency_cd=USGS&format=rdb'.format(gage)
    df = pd.read_csv(url, skiprows=64, sep = '\t')
    df.drop(0, axis=0, inplace=True)
    y = df['peak_va'].astype(float)
    x = df['peak_dt']
    x = pd.to_datetime(x, format= '%Y-%m-%d',  errors='coerce')
    y.index = x
    return y     

def MergeDatsets(tsvs):
    f = tsvs[0]           # Used the first data file to initialize the dataframe
    flow_type = f[-9:].split('.tsv')[0]
    print(os.path.basename(f))

    if flow_type == '60_dv':
        df = read_usgs_60_dv(f)
        
    elif flow_type == '60_iv':
        df = read_usgs_60_iv(f)

    elif flow_type == '65_iv':
        df = read_usgs_65_iv(f)

    else:
        print("Error")

    for f in tsvs[1:]:
        print(os.path.basename(f))
        data = f[-9:].split('.tsv')[0]
        gage =os.path.basename(f).split('.')[0]
        
        if data == '60_dv':
            df1 = read_usgs_60_dv(f)
            df = pd.merge(df, df1, how='outer', left_index=True, right_index=True)
            
        elif data == '60_iv':
            df1 = read_usgs_60_iv(f)
            df = pd.merge(df, df1, how='outer', left_index=True, right_index=True)
            
        elif data == '65_iv':
            df1 = read_usgs_65_iv(f)
            df = pd.merge(df, df1, how='outer', left_index=True, right_index=True)
            
        else:
            continue

    return df

def PlotPeaks(df, gage):
    ax1 = plt.subplot(111)
    df.sort_values(axis=0, ascending=False, inplace=True)
    plt.scatter(df.index , df, color = 'b', marker = 'o', facecolors='none', s=20)
    plt.title(gage.split('_')[0] + ' Peak Flows')
    plt.grid()
'''
def Plot_iv_Peak(df, peaks, rank, gage, spline=False, days=10):
    epsilon = pd.Timedelta('{} day'.format(days))
    peak = peaks.index[rank]
    start = peak - epsilon 
    stop = peak + epsilon
    start = start.strftime(format = '%Y-%m-%d')
    stop = stop.strftime(format = '%Y-%m-%d')
    data = df[gage][start:stop].interpolate()
    data.plot(figsize=(10,2), grid='on', title = gage, marker='o', markersize= 1.5)

def Plot_dv_Peak(df, peaks, rank, gage, spline=False, days = 10):
    epsilon = pd.Timedelta('{} day'.format(days))
    peak = peaks.index[rank]
    start = peak - epsilon 
    stop = peak + epsilon
    start = start.strftime(format = '%Y-%m-%d')
    stop = stop.strftime(format = '%Y-%m-%d')
    if spline:
        data = df[gage][start:stop].interpolate(method='spline', order=2)
    else:
        data = df[gage][start:stop].interpolate()
    data.plot(figsize=(10,2), grid='on', title = gage, marker='o', markersize= 1.5) 
'''
def PlotPeak(df, peaks, rank, gage, spline=False, days=10):
    epsilon = pd.Timedelta('{} day'.format(days))
    peak = peaks.index[rank]
    start = peak - epsilon 
    stop = peak + epsilon
    start = start.strftime(format = '%Y-%m-%d')
    stop = stop.strftime(format = '%Y-%m-%d')
    if spline:
        data = df[gage][start:stop].interpolate(method='spline', order=2)
    else:
        data = df[gage][start:stop].interpolate()
    data.plot(figsize=(10,2), grid='on', title = gage, marker='o', markersize= 1.5) 

def Q_to_Stage(usgs_gage_id):
    url = r'https://waterdata.usgs.gov/nwisweb/get_ratings?site_no={}&file_type=exsa'.format(usgs_gage_id)
    cols = ['INDEP', 'SHIFT','DEP','STOR']
    usgs_data = pd.read_csv(url, skiprows = 38, sep = '\t', names = cols)

    poly_order = 3
    xs = np.array(usgs_data['DEP'])
    ys = np.array(usgs_data['INDEP'])

    coefs    = np.polyfit(xs,ys,poly_order)
    polynomial = np.poly1d(coefs) 
    
    return polynomial

def GetPKFQ(gage):
    url = 'https://nwis.waterdata.usgs.gov/ny/nwis/peak?site_no={}&agency_cd=USGS&format=hn2'.format(gage)
    pkf = pd.read_csv(url)
    pkf.to_csv('return_periods\\{}.pkf'.format(gage), sep='\t', index=False)
    print('{} Data Saved in return_periods'.format(gage))
    
def GotoUSGS(state):
    url = 'https://waterdata.usgs.gov/nwis/uv?referred_module=sw&state_cd={}&site_tp_cd=OC&site_tp_cd=OC-CO&site_tp_cd=ES&site_tp_cd=LK&site_tp_cd=ST&site_tp_cd=ST-CA&site_tp_cd=ST-DCH&site_tp_cd=ST-TS&format=station_list'.format(state)
    print("\nCLICK HERE FOR USGS GAGES: \n", url)
    print("\nCLICK HERE FOR MAP: \n", 'https://maps.waterdata.usgs.gov/mapper/index.html')


def read_usgs_60_dv(f):
    gage = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, sep='\t')
    df = df.set_index(pd.DatetimeIndex(df['Date'])+ pd.Timedelta('12 hours'))
    df.rename(columns = {'Flow':gage}, inplace=True)
    df.drop(labels=['Date', 'agency_cd', 'site_no', 'Flow_cd'], axis = 1, inplace=True)
    return df

def read_usgs_60_iv(f):
    gage = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, sep='\t')
    df = df.set_index(pd.DatetimeIndex(df['dateTime']))
    df.rename(columns = {'Flow_Inst':gage}, inplace=True)
    df.drop(axis=1,labels =['dateTime','agency_cd','site_no', 'Flow_Inst_cd', 'tz_cd'], inplace=True)
    return df

def read_usgs_65_iv(f):
    gage = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, sep='\t')
    df = df.set_index(pd.DatetimeIndex(df['dateTime']))
    df.rename(columns = {'GH_Inst':gage}, inplace=True)
    df.drop(axis=1,labels =['dateTime','agency_cd','site_no', 'GH_Inst_cd', 'tz_cd'], inplace=True)
    return df



def CompareVolumes(instantaneous, daily, units = 60):
    '''
    Compare (unstretched) daily flows with instantaneous flows 
    '''
    # units for y values in cfs shuld be seconds, if resampled to 1 minute, value should be 60 
    ins = pd.DataFrame(instantaneous.copy())
    ins*=units

    dlymn = pd.DataFrame(daily.copy())
    dlymn*=units

    inst_volume, daily_volume = trapz(np.array(ins), x=None, dx=1.0, axis=0), trapz(np.array(dlymn), x=None, dx=1.0, axis=0)
    print('Volume from Instantaneous Observations = \t{}'.format(int(inst_volume)))
    print('Volume from Daily Mean Observations = \t\t{}'.format(int(daily_volume)))
    print('\nUsing Daily means yields = {} more Cubic Feet of Water'.format(int(daily_volume-inst_volume)))
    print('(Daily means results in a difference of volume of ~ {}% )'.format(float(100*(daily_volume-inst_volume)/inst_volume)))


def CompareVolumes_stretched(instantaneous, daily, units = 60):
    '''
    Compare  stretched daily flows 
    '''
    # units for y values in cfs shuld be seconds, if resampled to 1 minute, value should be 60 
    ins = pd.DataFrame(instantaneous.copy())
    ins*=units

    dlymn = pd.DataFrame(daily.copy())
    dlymn*=units

    inst_volume, daily_volume = trapz(np.array(ins), x=None, dx=1.0, axis=0), trapz(np.array(dlymn), x=None, dx=1.0, axis=0)
    print('Volume from Instantaneous Observations = \t{}'.format(int(inst_volume)))
    print('Volume from Stretched Daily Mean Observations = \t\t{}'.format(int(daily_volume)))
    print('\nUsing Stretched Daily means yields = {} more Cubic Feet of Water'.format(int(daily_volume-inst_volume)))
    print('(Streteched Daily means results in a difference of volume of ~ {}% )'.format(float(100*(daily_volume-inst_volume)/inst_volume)))


def PlotCumIntegral(instantaneous, daily, units = 60):
    #--Plot Cumlative flow (integrated)
    f, ax = plt.subplots(figsize=(10,2))
    ax.set_title('Cumulative Trapezoid')

    x = np.array(daily)*60
    ct = cumtrapz(x)
    ax.plot(ct, label = 'daily')

    x = np.array(instantaneous)*60
    ct = cumtrapz(x)
    ax.plot(ct, label = 'inst')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.grid()