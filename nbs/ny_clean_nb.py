
import matplotlib.pyplot as plt
from importlib import reload
import utils; reload(utils)
from utils import *
import os
from glob import glob
from pathlib import Path


def initialize(plot=False):
    #--Read in gage data
    data_dir = os.path.join(str(Path(os.getcwd()).parents[0]),'sample_data')
    tsvs = glob(os.path.join(data_dir, '015*'))
    metadata = os.path.join(data_dir ,'gage_metadata.tsv')
    df = MergeDatsets(tsvs) 
    return df, data_dir


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

    return daily_peak, inst_peak
    
    print('Maximum Daily Flow = {} cfs'.format(round(daily_peak.max(),0)))
    print('Maximum Inst Flow = {} cfs'.format(round(inst_peak.max(),0)))


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
    f, (ax) = plt.subplots()
    ax.plot(hydro.index, stage, label = 'Flood Stage (ft)')
    ax.axhline(breach_stage, color = 'red', label = 'Breach Height')
    ax.legend(loc='best',  fontsize='small')
    f.set_size_inches(6,4)
    ax.grid()
    f.autofmt_xdate()
    stage = pd.DataFrame( stage, index = hydro.index)
    stage.rename(columns = {0: 'stage'}, inplace=True)
    return stage