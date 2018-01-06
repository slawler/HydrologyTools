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
from pathlib import Path
from scipy.integrate import trapz, cumtrapz, simps
import numpy as np
import os
from glob import glob
from IPython.display import Markdown, display

#-----------------------------------------------------------------------------#
#--Notebook Functions
#-----------------------------------------------------------------------------#

def printbold(string):
    # Jupyter Notebook function to print bold
    display(Markdown('**'+string+'**'))


#-----------------------------------------------------------------------------#
#--Data Initialization & Manipulation Functions
#-----------------------------------------------------------------------------#

def initialize_test_case(plot=False):
    #--Read in gage data for the Michigan Test Case
    data_dir = os.path.join(str(Path(os.getcwd()).parents[0]),'sample_data')
    tsvs = glob(os.path.join(data_dir, '04119000*'))   
    printbold("Read in USGS Gage Records")    
    df = MergeDatsets(tsvs) 
    return df, data_dir

def initialize(plot=False):
    #--Read in gage data for this region (gages = 015*)
    data_dir = os.path.join(str(Path(os.getcwd()).parents[0]),'sample_data')
    tsvs = glob(os.path.join(data_dir, '015*'))
    metadata = os.path.join(data_dir ,'gage_metadata.tsv')
    df = MergeDatsets(tsvs) 
    return df, data_dir

def MergeDatsets(tsvs):
    #--Read in all tsvs and create a single dataframe (pass a list of tsvs)
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

def read_usgs_60_dv(f):
    #--Read in to dataframe, and reformat daily flow records
    gage = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, sep='\t')
    df = df.set_index(pd.DatetimeIndex(df['Date'])+ pd.Timedelta('12 hours'))
    df.rename(columns = {'Flow':gage}, inplace=True)
    df.drop(labels=['Date', 'agency_cd', 'site_no', 'Flow_cd'], axis = 1, inplace=True)
    return df

def read_usgs_60_iv(f):
    #--Read in to dataframe, and reformat instantaneous flow records
    gage = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, sep='\t')
    df = df.set_index(pd.DatetimeIndex(df['dateTime']))
    df.rename(columns = {'Flow_Inst':gage}, inplace=True)
    df.drop(axis=1,labels =['dateTime','agency_cd','site_no', 'Flow_Inst_cd', 'tz_cd'], inplace=True)
    return df

def read_usgs_65_iv(f):
    #--Read in to dataframe, and reformat instantaneous stage records
    gage = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, sep='\t')
    df = df.set_index(pd.DatetimeIndex(df['dateTime']))
    df.rename(columns = {'GH_Inst':gage}, inplace=True)
    df.drop(axis=1,labels =['dateTime','agency_cd','site_no', 'GH_Inst_cd', 'tz_cd'], inplace=True)
    return df


#-----------------------------------------------------------------------------#
#--USGS API-like Functions
#-----------------------------------------------------------------------------#

def Get_USGS_Peaks(gage):  
    #--Get peaks for a gage, may be sensitive to skiprows value. (check if superceded by GetPKFQ below)
    url = 'https://nwis.waterdata.usgs.gov/ny/nwis/peak?site_no={}&agency_cd=USGS&format=rdb'.format(gage)
    df = pd.read_csv(url, skiprows=64, sep = '\t')
    df.drop(0, axis=0, inplace=True)
    y = df['peak_va'].astype(float)
    x = df['peak_dt']
    x = pd.to_datetime(x, format= '%Y-%m-%d',  errors='coerce')
    y.index = x
    return y     

def GetPKFQ(gage):
    #--Need to test this function (see Get_USGS_Peaks)
    url = 'https://nwis.waterdata.usgs.gov/ny/nwis/peak?site_no={}&agency_cd=USGS&format=hn2'.format(gage)
    pkf = pd.read_csv(url)
    pkf.to_csv('return_periods\\{}.pkf'.format(gage), sep='\t', index=False)
    print('{} Data Saved in return_periods'.format(gage))
    
def GotoUSGS(state):
    #--Pass state initials (e.g. 'NY'), function opens gage catalogue in browser
    #--URL stirng broken into parts for ease of reading in text editor
    part1 = 'https://waterdata.usgs.gov/nwis/uv?referred_module=sw&state_cd={}'.format(state)
    part2 = '&site_tp_cd=OC&site_tp_cd=OC-CO&site_tp_cd=ES&site_tp_cd=LK&site_tp_cd=ST&'
    part3 = 'site_tp_cd=ST-CA&site_tp_cd=ST-DCH&site_tp_cd=ST-TS&format=station_list'
    url = part1 + part2 + part3
    print("\nCLICK HERE FOR USGS GAGES: \n", url)
    print("\nCLICK HERE FOR MAP: \n", 'https://maps.waterdata.usgs.gov/mapper/index.html')


#-----------------------------------------------------------------------------#
#--Math Functions
#-----------------------------------------------------------------------------#

def Q_to_Stage(usgs_gage_id, poly_order = 3):
    #--Use USGS rating curve to convert flow to stage (where available)
    url = r'https://waterdata.usgs.gov/nwisweb/get_ratings?site_no={}&file_type=exsa'.format(usgs_gage_id)
    cols = ['INDEP', 'SHIFT','DEP','STOR']
    usgs_data = pd.read_csv(url, skiprows = 38, sep = '\t', names = cols)
    
    xs = np.array(usgs_data['DEP'])
    ys = np.array(usgs_data['INDEP'])

    coefs    = np.polyfit(xs,ys,poly_order)
    polynomial = np.poly1d(coefs) 
    
    return polynomial

def IntegrateHydrograph(timeseries, da_sqft, method = 'trapezoid'):
    from scipy import integrate
    #--Integrate the volume under a hydrograph
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


def GetBreachFlow(base_hydrograph, site, rasdata, station, breach_point, breach_elev, data_dir, date_int=6, poly_order = 3):
    #--Calls ComputeWeirFlow to calculate weir flow & plot (would benefit from partitioning into simpler functions)
    df = GetRasData(rasdata, station)

    # Plot a Rating Curve using Stage & Flow data
    StageDischargePlot(df, figsize=(5,3))
    polyfit  = np.polyfit(df['flow'],df['stage'],poly_order)
    polyfitline = np.poly1d(polyfit) 
    hydro =pd.DataFrame(base_hydrograph)
    colname = hydro.columns
    hydro.rename(columns ={colname[0]:'Base Hydrograph'}, inplace=True)
    stage = plotcomp(hydro, polyfitline,breach_elev)
    df_weir = ComputeWeirFlow(stage, breach_elev, date_int) 
    output_csv = os.path.join(data_dir, '{}_BreachData_{}_location_{}.tsv'.format(site, station, breach_point))
    df_weir.to_csv(output_csv, sep = '\t')
    
    printbold('\nInflow Data for Breach Location: ')
    print(output_csv)  

def ComputeWeirFlow(df, breach_height, date_int, weir_coeff=2.0, breach_length=250, plot=True):
    #--Compute & plot weir flow at given HEC-RAS cross section, return a plot-series 
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

#-----------------------------------------------------------------------------#
#--HEC-RAS UTILS
#-----------------------------------------------------------------------------#

def GetRasData(hdf_plan_file, station):
    #--Read HEC-RAS hdf plan file, return dataframe with flow & calculated water surface elevation at a given station
    with h5py.File(hdf_plan_file ,'r') as hf:
        #--Paths should be the same in any hdf plan file created by HEC-RAS up to 5.0.3
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
    #--Read in unsteady flow data from HEC-RAS plan to dataframe
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

#-----------------------------------------------------------------------------#
#--Plotting Functions
#-----------------------------------------------------------------------------#

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

def PlotPeaks(df, gage):
    ax1 = plt.subplot(111)
    df.sort_values(axis=0, ascending=False, inplace=True)
    plt.scatter(df.index , df, color = 'b', marker = 'o', facecolors='none', s=20)
    plt.title(gage.split('_')[0] + ' Peak Flows')
    plt.grid()

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

#-----------------------------------------------------------------------------#
#--REGION SPECIFIC FUNCTIONS (copied from ny_clean_nb.py, version control cleanup) 
#
#--Functions created on the fly with a specific set of variables (i.e. gage, flow type, date, etc.)
#--TODO: Clean up functions so that they can be used universally.
#-----------------------------------------------------------------------------#

def Daily_vs_Instant(df):
    print('\n')
    printbold("Interpolated daily flows result in a lower peak flow, but a greater total volume of water:")

    # Daily Plot
    gage = '01509000_00060_dv'
    start, stop = '2005-03-15','2005-04-10'
    adj_start, adj_stop = '2005-03-29','2005-04-10'

    daily_peak= df[gage][start:stop].copy()
    daily_peak= daily_peak[adj_start:adj_stop].resample('60S').asfreq()
    daily_peak= daily_peak.interpolate(method = 'spline', order = 2)
    daily_peak = daily_peak['2005-04-02':'2005-04-05 12:00']
    daily_peak.name = 'Daily Flow'
    daily_peak.plot(figsize=(10,2), grid='on')
    
    # Instantaneous Plot
    gage = '01509000_00060_iv'
    inst_peak= df[gage][start:stop].copy()
    inst_peak= inst_peak[adj_start:adj_stop].resample('60S').asfreq()

    inst_peak= inst_peak.interpolate(method = 'spline', order = 3)
    inst_peak = inst_peak['2005-04-02':'2005-04-05 12:00']
    inst_peak.name = 'Instantaneous Flow'
    inst_peak.plot(figsize=(10,2), grid='on')

    plt.legend()
    plt.title('USGS Gage {}'.format(gage.split('_')[0]))
    print('Maximum Daily Flow = {} cfs'.format(round(daily_peak.max(),0)))
    print('Maximum Inst Flow = {} cfs'.format(round(inst_peak.max(),0)))

    return daily_peak, inst_peak
    



def Stretched_Daily_vs_Instant(df):
    print('\n')
    printbold('Results from Stretched Daily flows compared with Instantaneous Records')
    # Daily Plot
    gage = '01509000_00060_dv'
    start, stop = '2005-03-15','2005-04-10'
    adj_start, adj_stop = '2005-03-29','2005-04-10'

    #--Replace the daily mean on the day of peak flow, with the peak from the instantaneous record
    daily_peak= df[gage][start:stop].copy()
    daily_peak['2005-04-03 12:00:00'] = 14200.0 
    daily_peak= daily_peak[adj_start:adj_stop].resample('60S').asfreq()
    daily_peak= daily_peak.interpolate(method = 'spline', order = 2)
    daily_peak = daily_peak['2005-04-02':'2005-04-05 12:00']
    daily_peak.name = 'Daily Flow, Stretched to Peak'
    daily_peak.plot(figsize=(10,2), grid='on')

    # Instantaneous Plot
    gage = '01509000_00060_iv'
    inst_peak= df[gage][start:stop].copy()
    inst_peak= inst_peak[adj_start:adj_stop].resample('60S').asfreq()
    inst_peak= inst_peak.interpolate(method = 'spline', order = 3)
    inst_peak = inst_peak['2005-04-02':'2005-04-05 12:00']
    inst_peak.name = 'InstantaneousFlow'
    inst_peak.plot(figsize=(10,2), grid='on')

    plt.title('{}'.format(gage.split('_')[0]))
    plt.legend()

    print('Maximum Daily Flow (stretched) = {} cfs'.format(round(daily_peak.max(),0)))
    print('Maximum Inst Flow = {} cfs'.format(round(inst_peak.max(),0)))
    return daily_peak, inst_peak


def Stretched_Daily_100yr(df, plot = False):
    print('\n')

    # Daily Plot
    gage = '01509000_00060_dv'
    start, stop = '2005-04-01','2005-04-07'

    #--Replace the daily mean on the day of peak flow, with the peak from the instantaneous record
    daily_peak= df[gage][start:stop].copy()
    daily_peak['2005-04-03 12:00:00'] = 20960.0 
    daily_peak= daily_peak.resample('15T').asfreq().interpolate(method = 'spline', order = 2)
    daily_peak.name = '1 Pct'
    daily_peak = daily_peak['2005-04-02':'2005-04-05 12:00']


    # Instantaneous Plot
    gage = '01509000_00060_iv'
    inst_peak= df[gage][start:stop].copy().interpolate(method = 'spline', order =2)
    inst_peak.name = 'Base Hydrograph'
    inst_peak = inst_peak['2005-04-02':'2005-04-05 12:00']
    if plot == True:
        daily_peak.plot(figsize=(8,6), grid='on',color = 'blue')
        inst_peak.plot(figsize=(8,6), grid='on', color = 'green')

        plt.title('{}'.format(gage.split('_')[0]))
        plt.legend(loc=0)
        printbold("Hydrograph properties:")
        print('1 Percent Peak Flow = {} cfs'.format(round(daily_peak.max(),0)))
        print('Maximum Inst Flow = {} cfs'.format(round(inst_peak.max(),0)))
    return daily_peak, inst_peak


#FINAL HYDROGRAPH:
def Hydrograph_100yr(df, peak):
    print('\n')
    printbold('Historic Hydrograph Adapted to 100yr Peak Flow')
    # Daily Plot
    gage = '01509000_00060_dv'
    start, stop = '2005-03-15','2005-04-10'
    adj_start, adj_stop = '2005-03-29','2005-04-10'

    #--Replace the daily mean on the day of peak flow, with the peak from the instantaneous record
    daily_peak= df[gage][start:stop].copy()
    daily_peak['2005-04-03 12:00:00'] = peak 
    daily_peak= daily_peak[adj_start:adj_stop].resample('60S').asfreq()
    daily_peak= daily_peak.interpolate(method = 'spline', order = 2)
    daily_peak = daily_peak[3600:-3600]
    daily_peak.name = 'Daily Values, Stretched to Inst. Peak'
    daily_peak.plot(figsize=(12,3), grid='on')

    plt.title('{}'.format(gage.split('_')[0]))
    plt.legend(location = 0)
    
    return pd.DataFrame(daily_peak)

def plotcomp(hydro, rawfitline, breach_stage=0):
    stage = rawfitline(hydro['Base Hydrograph'])
    f, (ax, ax1) = plt.subplots(1,2, figsize = (10,3))
    ax.plot(hydro.index, hydro['Base Hydrograph'], label = '1 Pct Flow (cfs)')
    ax.legend(loc='best',  fontsize='small')
    ax.grid()

    ax1.plot(hydro.index, stage, label = 'Flood Stage (ft)')
    ax1.axhline(breach_stage, color = 'red', label = 'Breach Height')
    ax1.legend(loc='best',  fontsize='small')
    ax1.grid()

    #--Format the Plots
    import matplotlib.dates as mdates
    ax.set_title('Flow Hydrograph \n (input from smoothed storm )')
    days = mdates.DayLocator()   # every year
    dateFmt = mdates.DateFormatter('%Y-%m-%d')
    ax.xaxis.set_major_formatter(dateFmt)
    ax.xaxis.set_major_locator(days)

    ax1.set_title('Stage Hydrograph \n(output using rating curve above)')
    ax1.xaxis.set_major_locator(days)
    ax1.xaxis.set_major_formatter(dateFmt)    
    f.autofmt_xdate()

    stage = pd.DataFrame( stage, index = hydro.index)
    stage.rename(columns = {0: 'stage'}, inplace=True)
    return stage

def init_base_hydro(gage_data, peak=14209.0, pct_1=20960.0):
    pct_1_peak, inst_peak  = Stretched_Daily_100yr(gage_data, plot = False)
    stretch_1pct = pct_1/peak 
    smooth_storm_resample = inst_peak.resample('30T').mean()
    smooth_storm_1pct = smooth_storm_resample*stretch_1pct
    #print('Factor: \n','\t1 Percent\t{}'.format(stretch_1pct))

    f, ax = plt.subplots()
    ax.plot(smooth_storm_1pct,color = 'green',label = 'Base Hydrograph from Gage, stretched to 1-Pct')
    ax.grid()
    ax.legend(loc='best')
    title = '{} - {}'.format(smooth_storm_1pct.index[0],smooth_storm_1pct.index[-1])
    ax.set_title(title)

    f.set_size_inches(6,4)
    f.autofmt_xdate()
    return smooth_storm_1pct

def smooth_base_hydro(smooth_storm_1pct):
    # 1. Make a copy of the unsmoothed hydrograph: 
    final_hydrograph  = smooth_storm_1pct.copy()

    #2. Replace the dates of the portions needing to be smoothed with NaNs
    final_hydrograph['2005-04-02 00:00':'2005-04-03 00:00'] = np.nan #

    #3. Set the first data point to the desired starting value 
    # Note: Typically the starting and ending values will be the same, for this example the min works
    final_hydrograph['2005-04-02 00:00':'2005-04-02 00:00'] = float(final_hydrograph.min())

    #4. Add more points as needed to smooth the curve
    # Note: This is an iterative process, using a (factor*min) appproach is used here for speed of iteration 
    final_hydrograph['2005-04-02 06:00':'2005-04-02 07:00'] = float(final_hydrograph.min()*1.05)
    final_hydrograph['2005-04-02 14:00':'2005-04-02 15:00'] = float(final_hydrograph.min()*1.25)

    #5. Smooth the iterated points
    smooth_storm = final_hydrograph.interpolate(how = 'polynomial', order=3) 

    f, ax = plt.subplots()
    ax.plot(smooth_storm,color = 'green',label = 'Smoothed Storm')
    ax.grid()
    ax.legend(loc= 'best',  fontsize='small')

    f.set_size_inches(6,4)
    f.autofmt_xdate()
    return smooth_storm
  