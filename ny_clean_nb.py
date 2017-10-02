
import matplotlib.pyplot as plt
from importlib import reload
import utils; reload(utils)
from utils import *
import os
from glob import glob


def initialize():
    #--Read in gage data
    data_dir =os.getcwd()
    region_name = "TIOUGHNIOGA_RIVER"
    output_dataset = os.path.join(data_dir, '{}.pkl'.format(region_name))
    tsvs = glob(os.path.join(data_dir, '*.tsv'))

    metadata = os.path.join(data_dir ,'gage_metadata.tsv')

    if metadata in tsvs:
        tsvs.remove(metadata)
    #print('\n')    
    printbold("Read in USGS Gage Records")    
    df = MergeDatsets(tsvs) 

    #df.head().style    
    gage = '01509520_00065_iv'
    df[gage]['2017-05-25':'2017-10'].plot(figsize=(10,2), grid='on', title = 'Stage (ft) for Period of Record \n USGS {}'.format(gage.split('_')[0]))
    return df


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
    daily_peak = daily_peak[3600:-3600]
    daily_peak.name = 'Daily Peak Flow, Interpolated'
    daily_peak.plot(figsize=(10,2), grid='on')
    
    # Instantaneous Plot
    gage = '01509000_00060_iv'
    inst_peak= df[gage][start:stop].copy()
    inst_peak= inst_peak[adj_start:adj_stop].resample('60S').asfreq()

    inst_peak= inst_peak.interpolate(method = 'spline', order = 3)
    inst_peak = inst_peak[3600:-3600]
    inst_peak.name = 'Instantaneous Peak Flow, Interpolated'
    inst_peak.plot(figsize=(10,2), grid='on')

    plt.legend()
    plt.title('USGS Gage {}'.format(gage.split('_')[0]))
    
    print('Maximum Daily Flow = {} cfs'.format(round(daily_peak.max(),0)))
    print('Maximum Inst Flow = {} cfs'.format(round(inst_peak.max(),0)))
    CompareVolumes(inst_peak, daily_peak)
    PlotCumIntegral(inst_peak, daily_peak)


def Stretched_Daily_vs_Instant(df, plot_cumulative=False):
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
    daily_peak = daily_peak[3600:-3600]
    daily_peak.name = 'Daily Values, Stretched to Inst. Peak'
    daily_peak.plot(figsize=(10,2), grid='on')

    # Instantaneous Plot
    gage = '01509000_00060_iv'
    inst_peak= df[gage][start:stop].copy()
    inst_peak= inst_peak[adj_start:adj_stop].resample('60S').asfreq()
    inst_peak= inst_peak.interpolate(method = 'spline', order = 3)
    inst_peak = inst_peak[3600:-3600]
    inst_peak.name = 'Instantaneous Peak Flow, Interpolated'
    inst_peak.plot(figsize=(10,2), grid='on')

    plt.title('{}'.format(gage.split('_')[0]))
    plt.legend()

    print('Maximum Daily Flow (stretched) = {} cfs'.format(round(daily_peak.max(),0)))
    print('Maximum Inst Flow = {} cfs'.format(round(inst_peak.max(),0)))
    CompareVolumes_stretched(inst_peak, daily_peak)
    if plot_cumulative:
        PlotCumIntegral(inst_peak, daily_peak)



def Stretched_Daily_100yr(df, plot_cumulative=False):
    print('\n')
    printbold("Hydrograph properties:")
    # Daily Plot
    gage = '01509000_00060_dv'
    start, stop = '2005-03-15','2005-04-10'
    adj_start, adj_stop = '2005-03-30','2005-04-09'

    #--Replace the daily mean on the day of peak flow, with the peak from the instantaneous record
    daily_peak= df[gage][start:stop].copy()
    daily_peak['2005-04-03 12:00:00'] = 20960.0 
    daily_peak= daily_peak[adj_start:adj_stop].resample('60S').asfreq()
    daily_peak= daily_peak.interpolate(method = 'spline', order = 2)
    daily_peak = daily_peak[3600:-3600]
    daily_peak.name = 'Daily Values, Stretched to 100yr'
    daily_peak.plot(figsize=(10,2), grid='on')

    # Instantaneous Plot
    gage = '01509000_00060_iv'
    inst_peak= df[gage][start:stop].copy()
    inst_peak= inst_peak[adj_start:adj_stop].resample('60S').asfreq()
    inst_peak= inst_peak.interpolate(method = 'spline', order = 3)
    inst_peak = inst_peak[3600:-3600]
    inst_peak.name = 'Instantaneous Peak Flow, Interpolated'
    inst_peak.plot(figsize=(10,2), grid='on')

    plt.title('{}'.format(gage.split('_')[0]))
    plt.legend()

    print('Maximum Daily Flow (stretched) = {} cfs'.format(round(daily_peak.max(),0)))
    print('Maximum Inst Flow = {} cfs'.format(round(inst_peak.max(),0)))
    CompareVolumes_stretched(inst_peak, daily_peak)
    if plot_cumulative:
        PlotCumIntegral(inst_peak, daily_peak)


def Stretched_Daily_100yr_2(df, plot_cumulative=False):
    print('\n')
    printbold("Hydrograph properties:")
    # Daily Plot
    gage = '01509000_00060_dv'
    start, stop = '1979-02-15','1979-03-15'
    adj_start, adj_stop = '1979-03-01','1979-03-12'

    #--Replace the daily mean on the day of peak flow, with the peak from the instantaneous record
    daily_peak= df[gage][start:stop].copy()
    daily_peak['1979-03-06 12:00:00'] = 20960.0 
    daily_peak= daily_peak[adj_start:adj_stop].resample('60S').asfreq()
    daily_peak= daily_peak.interpolate(method = 'spline', order = 3)
    daily_peak = daily_peak[3600:-3600]
    daily_peak.name = 'Daily Values, Stretched to 100yr'
    daily_peak.plot(figsize=(10,2), grid='on')

    # Instantaneous Plot
    inst_peak= df[gage][start:stop].copy()
    inst_peak= inst_peak[adj_start:adj_stop].resample('60S').asfreq()
    inst_peak= inst_peak.interpolate(method = 'spline', order = 3)
    inst_peak = inst_peak[3600:-3600]
    inst_peak.name = 'Instantaneous Peak Flow, Interpolated'
    inst_peak.plot(figsize=(10,2), grid='on')

    plt.title('{}'.format(gage.split('_')[0]))
    plt.legend()

    print('Maximum Daily Flow (stretched) = {} cfs'.format(round(daily_peak.max(),0)))
    print('Maximum Inst Flow = {} cfs'.format(round(inst_peak.max(),0)))
    CompareVolumes_stretched(inst_peak, daily_peak)
    if plot_cumulative:
        PlotCumIntegral(inst_peak, daily_peak)
