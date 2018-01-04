# -*- coding: utf-8 -*-
"""
@author: slawler
"""
import pandas as pd
import h5py
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
import os
from glob import glob
from IPython.display import Markdown, display

def initialize_test_case(plot=False):
    #--Read in gage data
    data_dir = os.path.join(str(Path(os.getcwd()).parents[0]),'sample_data')
    tsvs = glob(os.path.join(data_dir, '04119000*'))
    #print('\n')    
    printbold("Read in USGS Gage Records")    
    df = MergeDatsets(tsvs) 
    return df, data_dir

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
def PlotMultiPeak(df, peaks, rank, gage, spline=False, days=10):
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
    plot = data.plot(figsize=(10,2), grid='on', title = gage, marker='o', markersize= 1.5)  
    plot.set_xlim(pd.datetime(1900, 1, 1) , pd.datetime(2020, 1, 1))
    return plot

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
    plot = data.plot(figsize=(10,2), grid='on', title = gage, marker='o', markersize= 1.5)  
    return plot

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

def IntegrateHydrograph(timeseries, da_sqft, method = 'trapezoid'):
    from scipy import integrate
    '''Flow should be given in units of cfs, index should be time series'''
    data = timeseries.resample('1S').asfreq().interpolate()
    da_acres = 640*da_sqft
    
    if method == 'simpson':
        volume_acre_inches = volume*0.00027548209085905
        inches = volume_acre_inches/da_acres
        volume = integrate.simps(np.array(data.values), x=None, dx=1, axis=0, even='avg')
    else:
        volume = trapz(np.array(data.values), x=None, dx=1.0, axis=0)
        volume_acre_inches = volume*0.00027548209085905
        inches = volume_acre_inches/da_acres
        
    return round(float(inches), 2)

'''
!NEED TO UPDATE FUNCTION USING IntegrateHydrograph Approach
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
'''

#--RAS UTILS
def GetRasData(hdf_plan_file, station):
    with h5py.File(hdf_plan_file ,'r') as hf:
        xs = hf.get('/Results/Steady/Output/Geometry Info/Cross Section Only')[:]
        xs = pd.DataFrame([x.astype(str).split(' ') for x in xs], columns=['River', 'Reach', 'Station'])

        flows = hf.get('/Results/Steady/Output/Output Blocks/Base Output/Steady Profiles/Cross Sections/Flow')[:]
        flows = pd.DataFrame(flows)

        stages = hf.get('/Results/Steady/Output/Output Blocks/Base Output/Steady Profiles/Cross Sections/Water Surface')[:]
        stages = pd.DataFrame(stages)
        
    idx = xs[xs['Station'] == str(station)].index[0]
    q = flows[idx]
    s = stages[idx]

    df = pd.DataFrame({'flow':q, 'stage':s})
    df.name = station   
    df.sort_values(by='flow', inplace=True)
    printbold('RAS Data for XS {}'.format(df.name))
    return df

def GetRasUnsteadyFlow(hdf_plan_file):
    with h5py.File(hdf_plan_file ,'r') as hf:
        qs = hf.get('/Event Conditions/Unsteady/Boundary Conditions/Flow Hydrographs/')
        df=pd.DataFrame()
        for q in qs:
            flow_name = q.split(': ')[-1]
            print(flow_name)
            ats =  qs[q].attrs
            items = ats.items
            for item in items():

                if item[0] == 'Start Date': 
                    start = item[1]
                    start = start.astype(str)
                    if start.split(' ')[-1][:4] == '2400': #HEC-RAS likes to end the day with 2400...confuses python!
                        s = start.replace('2400', '0000')
                    start = datetime.strptime(s, '%d%b%Y %H%M')
                    #print('start', s)

                elif item[0] == 'End Date' :
                    end= item[1]
                    end = end.astype(str)
                    if end.split(' ')[-1][:4] == '2400':
                        e = end.replace('2400', '0000')
                    end = datetime.strptime(e, '%d%b%Y %H%M')
                    #print('end',e )
                else:
                    continue

            flow = qs[q][:]
            flow = flow[:,1]
            steps = len(flow) -1
            freq = (end-start)/steps
            idx = pd.date_range(start,end, freq=freq)
            df[flow_name] = flow
        df = df.set_index(pd.DatetimeIndex(idx))
        return df


def StageDischargePlot(df, figsize = (2,3)):
    x, y = df['flow'], df['stage']
    f, (ax) = plt.subplots()
    ax.plot(x,y)
    ax.grid()
    ax.set_ylabel('Stage (ft)')
    ax.set_xlabel('Flow (cfs)')
    ax.set_xlim(0, x.max()*1.1)
    ax.set_title('Rating Curve for XS {}'.format(df.name))
    f.set_size_inches(figsize)
    f.autofmt_xdate()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)



def ComputeWeirFlow(df, breach_height, date_int, weir_coeff=2.0, breach_length=250, plot=True):
    import matplotlib.dates as mdates
    min_stage = float(df.stage.min())
    df['head'] = df['stage']-breach_height
    df = df.query('head > 0 or head == 0')
    df['weir_flow'] = weir_coeff*breach_length*df['head']**(2/3)
    if plot: 
        f, (ax1, ax2, ax3) = plt.subplots(1,3, figsize = (12,5))
        
        ax1.plot(df['stage'],color = 'green',label = 'River Stage (ft) at Breach Location')
        ax1.grid()
        ax1.legend(loc= 'best',  fontsize='x-small')
        
        ax2.set_title('Normalized Hydrographs based on Breach Elevation \n')
        ax2.plot(df['head'],color = 'green',label = 'Head (ft) at Breach Location')
        ax2.grid()
        ax2.legend(loc= 'best',  fontsize='x-small')
        
        ax3.plot(df['weir_flow'],color = 'blue',label = 'Weir Flow (cfs)')
        ax3.grid()
        ax3.legend(loc= 'best',  fontsize='small')

        days = mdates.HourLocator(interval = date_int)   # every year
        dateFmt = mdates.DateFormatter('%Y-%m-%d :%H:%M')
        ax1.xaxis.set_major_formatter(dateFmt)
        ax1.xaxis.set_major_locator(days)
        ax2.xaxis.set_major_formatter(dateFmt)
        ax2.xaxis.set_major_locator(days)
        ax3.xaxis.set_major_formatter(dateFmt)
        ax3.xaxis.set_major_locator(days)

        f.autofmt_xdate()
        
    
    return pd.DataFrame(df)

def GetRasWSE(hdf_plan_file, profile):
    with h5py.File(hdf_plan_file ,'r') as hf:
        xs = hf.get('/Results/Steady/Output/Geometry Info/Cross Section Only')[:]
        xs = pd.DataFrame([x.astype(str).split(' ') for x in xs], columns=['River', 'Reach', 'Station'])   
        
        profiles = hf.get('/Results/Steady/Output/Output Blocks/Base Output/Steady Profiles/Profile Names')[:]
        profiles = pd.DataFrame({'p': profiles[:]})
        profiles['p'] = profiles['p'].str.decode("utf-8")
        idx = profiles[profiles['p'] == profile].index[0]

        stages = hf.get('/Results/Steady/Output/Output Blocks/Base Output/Steady Profiles/Cross Sections/Water Surface')[:]
        stages = pd.DataFrame(stages)
        stages = pd.DataFrame({profile: stages.iloc[0]})
        
    df = pd.merge(xs, stages, left_index=True, right_index=True, how = 'inner')
    return df