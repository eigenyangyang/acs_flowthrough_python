#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Underway flow-through ac-s data analysis.

Before using this script, please
1) prepare the following data:
- raw ac-s ASCII data exported from Wetview or Compass software;
- ac-s device file (*.dev);
- coincident tempreture and salinity data file (*.txt):
first column: date time in the format of "yyyy-mm-dd HH:MM:SS"; 
second column: temperature in degrees Celsius;
third column: salinity;
Separators between columns are commas(",");
- tempreture and salinity correction coefficients (e.g. Sullivan et al., 2006) 
file (*.xls).

2) install:
- Python 3.7x or higher version
- Essential packages: numpy, matplotlib, pandas, scipy.


Analysis method detailed in:
1) Liu et al. (2018) Underway spectrophotometry in the Fram Strait (European 
Arctic Ocean): a highly resolved chlorophyll a data source for complementing 
satellite ocean color. Optics Express, 26(14), A678-A696.
2) Liu et al. (2019) Retrieval of Phytoplankton Pigments from Underway
Spectrophotometry in the Fram Strait. Remote Sensing, 11(3), 318.

@author: Yangyang Liu (yangyang.liu@awi.de), June 2019.
'''

#-----------------------------------------------------------------------------
# Setup
#-----------------------------------------------------------------------------

import glob, os, sys, math
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as inp
from pandas import Series, DataFrame
from dateutil import parser #Date parser
import datetime
import scipy.optimize


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------     
def acs_rd_wetview(filename):
    
    '''
    This function extracts ac-s (Seabird Inc.) data from Wetview software 
    output (.dat files).
    
    Input: 
    filename - absolute path of the Wetview output raw ac-s file.
    
    Output: 
    1) Extracted ac-s data file saved as "*extracted.txt", in which the columns 
    are datetime, a and c, respectively.
    2) Extracted wavelength data saved as "wavelength.txt".
    '''
    
    foo = filename.split('/')
    print(f'Importing {foo[-1]}')
      
    try:
        
        with open(filename, 'r') as f:
            date_str = f.readline()
            date_str = date_str.split('\t') 
            date_str = date_str[1]+'\t'+date_str[2]
            date = parser.parse(date_str)
            date = pd.Series(date)
            print ('File recording starting time: ' + date_str)  
            for num, line in enumerate(f,1):
                if 'acquisition binsize' in line:
                    break   
            line = f.readlines()      
        line = line[0].split('\t')  
        
        wlc=list()
        wla=list()
        for wavelength in line:
            if wavelength.startswith(('C', 'c')):
                wl = float(wavelength[1:])
                wlc.append(wl)
  
            if wavelength.startswith(('A', 'a')):
                wl = float(wavelength[1:])
                wla.append(wl)
                
        if not len(wlc)==len(wla):
            print('Unexpected number of wavelength.\n')
        
        wla = pd.Series(wla)
        if not os.path.isfile('wavelength.txt'):
            wla.to_csv('wavelength.txt', sep='\t', index=False, 
                                header=False, encoding='utf-8')

        data = pd.read_csv(filename, header=None, skiprows=num+2, sep='\t') 
        data_dt = data.iloc[:,0]
        data_c = data.iloc[:,1:len(wla)+1]
        data_a = data.iloc[:,len(wla)+1:2*len(wla)+1]
        func = inp.interp1d(wlc, data_c, 'linear', fill_value='extrapolate')
        data_c = func(wla) #Interpolate data_c on wavelength_a scale
        data_c = DataFrame(data_c)
        
        #Check if there is a jump in timestmap of wetview
        i = data_dt.idxmax()
        if i != len(data_dt)-1:
            delta = np.median(data_dt.iloc[1:i+1] - data_dt.iloc[0:i])
            data_dt.iloc[i:-1] = data_dt.iloc[i:-1] + data_dt.iloc[i-1
                        ] - data_dt.iloc[i] + delta
        
    
        data_time = pd.Series()
        time_step = Series.tolist(data_dt)
        for j in range(len(time_step)):
            time_counter = date + datetime.timedelta(milliseconds=time_step[j])
            data_time = data_time.append(time_counter)

        data_time.index = range(len(data_time)) 
    
        data_t = data_time.apply(lambda x: x.strftime('%Y-%m-%d %H:%M:%S'))
    
        #Write output
        output = np.hstack([data_t.values[:,np.newaxis],
                            data_a.values, data_c.values])
            
        pd.DataFrame(data=output, columns=range(2*len(wla)+1)).to_csv(
                filename.replace('.dat','_extracted.txt'), index=False, 
                header=False, encoding='utf-8')
             
        print(f'Extraction of {foo[-1]} Finished')
    
    except FileNotFoundError:
        print(f'Error: Unable to open file: {filename}')


#-----------------------------------------------------------------------------    
#-----------------------------------------------------------------------------     
def acs_rd_compass(filename):
    
    '''
    This function extracts ac-s (Seabird Inc.) data from Compass software 
    output (.dat files).
    
    Input: 
    filename - absolute path of the Compass output raw ac-s file.
    
    Output: 
    1) Extracted ac-s data file saved as "*extracted.txt", in which the columns 
    are datetime, a and c, respectively.
    2) Extracted wavelength data saved as "wavelength.txt".
    '''
    
    foo = filename.split('/')
    print(f'Importing {foo[-1]}')
      
    try:
        
        with open(filename, 'r') as f:
            for num, line in enumerate(f,1):
                if 'Time(ms)' in line:
                    break      
        line = line.split('\t')       
    
        wlc=list()
        wla=list()
        for wavelength in line:
            if wavelength.startswith(('c', 'C')):
                wl = float(wavelength[1:])
                wlc.append(wl)
  
            if wavelength.startswith(('a', 'A')):
                wl = float(wavelength[1:])
                wla.append(wl)
                
        if not len(wlc)==len(wla):
            print('Unexpected number of wavelength.\n')
        
        wla = pd.Series(wla)
        if not os.path.isfile('wavelength.txt'):
            wla.to_csv('wavelength.txt', sep='\t', index=False, 
                                header=False, encoding='utf-8')

        data = pd.read_csv(filename, header=None, skiprows=num, sep='\t') 
        data_dt = data.iloc[:,0]
        data_c = data.iloc[:,1:len(wla)+1]
        data_a = data.iloc[:,len(wla)+1:2*len(wla)+1]
        func = inp.interp1d(wlc, data_c, 'linear', fill_value='extrapolate')
        data_c = func(wla) #Interpolate data_c on wavelength_a scale
        data_c = DataFrame(data_c)
    
        #Check if there is a jump in timestmap of compass
        i = data_dt.idxmax()
        if i != len(data_dt)-1:
            delta = np.median(data_dt.iloc[1:i+1] - data_dt.iloc[0:i])
            data_dt.iloc[i:-1] = data_dt.iloc[i:-1] + data_dt.iloc[i-1
                        ] - data_dt.iloc[i] + delta
        
        #Update date & time    
        date_str = foo[-1].split('_')[1][:-4]
        date = parser.parse(date_str)
        date = pd.Series(date)
        print ('File recording end time: ' + date_str)
    
        data_time = pd.Series()
        time_step = data_dt-data_dt.iloc[0]
        time_step = Series.tolist(time_step)
        for j in range(len(time_step)):
            time_counter = date + datetime.timedelta(milliseconds=time_step[j])
            data_time = data_time.append(time_counter)
    
        #File time stamp is done at the end of the recording, so we need to 
        #substract the length of the recording to the timestamp
        data_time = data_time - (data_time.iloc[-1]-data_time.iloc[0]);
        data_time.index = range(len(data_time)) 
    
        data_t = data_time.apply(lambda x: x.strftime('%Y-%m-%d %H:%M:%S'))
    
        #Write output
        output = np.hstack([data_t.values[:,np.newaxis],
                            data_a.values, data_c.values])
            
        pd.DataFrame(data=output, columns=range(2*len(wla)+1)).to_csv(
                filename.replace('.dat','_extracted.txt'), index=False, 
                header=False, encoding='utf-8')
             
        print(f'Extraction of {foo[-1]} Finished')
    
    except FileNotFoundError:
        print(f'Error: Unable to open file: {filename}')
    
    
#----------------------------------------------------------------------------- 
#----------------------------------------------------------------------------- 
def rsd_median(df):
    
    '''
    This function calculates relative median standard deviation (i.e. relative 
    median coefficient of variation, median standard deviation/median) of a 
    Pandas DataFrame for several rows.
    
    Input: df - Pandas DataFrame.
    
    Output: 
    relative median coefficient of variation on row axis.
    '''
    
    row = df.shape[0]
    df_median = np.nanmedian(df,axis=0)
    diff_squre=(df-df_median)**2 
    sumofsquare = np.nansum(diff_squre, axis=0)
    
    if row>1:
        sd_median = sumofsquare/(row-1)
    else:
        sd_median = sumofsquare/row

    rsd = np.sqrt(sd_median)/abs(df_median)
    
    return rsd
    

#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------- 
def acs_rm_spikes(filename, a_Thold, c_Thold, a_cvThold=None, c_cvThold=None):
    
    '''
    This function removes spikes of ac-s data. It combines 2 ways:
    1) ac data measurements whose data at 440nm are greater than the thresholds 
    a(c)_Thold are removed.
    2) ac data measured within 2 seconds whose relative median coefficient of 
    variations are greater than a(c)_cvThold are removed. Also, the 10 seconds 
    meaurements centered on these outliers are removed. The choices of the 2 
    thresholds relies on the understanding of the data. If the users are not 
    sure which values to set to them, then leave them as default, i.e. 
    a_cvThold=None, c_cvThold=None.

    Input: 
    filename - absolute path of the extracted raw ac-s file named as "*extracted.txt".
    a(c)_cvThold - threshold for the median coefficient of variation of 
    a(c) measurements in 2 seconds (8 measurements). e.g. 5.
    a(c)_Thold - threshold for a(c) data at 440nm, [m^-1]. e.g. 2.
    
    Output: despiked ac-s data saved as "*despiked.txt", in which the columns 
    are datetime, a and c, respectively.
    '''
    
    foo = filename.split('/')
    print(f'Removing spikes of {foo[-1]}')
    
    nmeas_1s = 4 #4 measurements per second.
    nmeas_10s = 10*4 #10 seconds' measurements
    
    wl = pd.read_csv('wavelength.txt', header=None, sep='\t')
    n_wl = len(wl)
    pos_NIR = np.where(wl>=700)[0]
    pos_440 = np.where(wl>=440)[0][0]

    try:
        
        #Read data file
        data = pd.read_csv(filename, header=None, sep=',')
        data_a = data.iloc[:,1:n_wl+1]
        data_c = data.iloc[:,n_wl+1:]
        data_t = data.iloc[:,0]
         
        pos = np.where((data_a.iloc[:,pos_440]>a_Thold) | 
               (data_c.iloc[:,pos_440]>c_Thold))[0]
        data_a.iloc[pos,:] = np.nan
        data_c.iloc[pos,:] = np.nan
         
        #Calculate median coefficient of variation for a and c    
        if a_cvThold is not None or c_cvThold is not None:
            for row in range(len(data_t)):     
                if row - nmeas_1s >= 0 and row + nmeas_1s <= len(data_t)-1:
                    nn = np.int_(range(row - nmeas_1s,row + nmeas_1s + 1))
                elif row - nmeas_1s < 0:
                    nn = np.int_(range(0,row + nmeas_1s+1))
                else:
                    nn = np.int_(range(row - nmeas_1s,len(data_t)))
                             
                a_cv = rsd_median(data_a.iloc[nn,range(pos_NIR[1])])
                c_cv = rsd_median(data_c.iloc[nn,range(pos_NIR[1])])
             
                if any(a_cv>a_cvThold) or any(c_cv>c_cvThold):
                
                    data_a.iloc[nn,:] = np.nan
                    data_c.iloc[nn,:] = np.nan
                
                    if row - nmeas_10s >= 0 and row + nmeas_10s <= len(data_t):
                        data_a.iloc[row - nmeas_10s : row + nmeas_10s,:] = np.nan
                        data_c.iloc[row - nmeas_10s : row + nmeas_10s,:] = np.nan
        
        #Write output
        output = np.hstack([data_t.values[:,np.newaxis],
                            data_a.values, data_c.values])
            
        pd.DataFrame(data=output, columns=range(2*n_wl+1)).to_csv(
                filename.replace('.txt','_despiked.txt'), index=False, 
                header=False, encoding='utf-8')
        
        print(f'Spikes removal of {foo[-1]} Finished!')     
     
        
    except FileNotFoundError:
        print(f'Error: Unable to open file: {filename}')


#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------- 
def acs_filemerge_oneminbin(dir_acs):
    
    '''
    This function merges all despiked ac data and bins them to 1-minute interval.
    
    Input:
    dir_acs - absolute path of the directory containing all despiked ac-s data 
    files "*despkied.txt".
    
    Output - merged and 1-min binned despiked ac-s data saved as 
    "acs_extracted_despiked_merged_1min.txt", in which the columns are datetime, 
    a and c, respectively.
    '''
    
    print ('Merging data files...')
    fnames=sorted(glob.glob(os.path.join(dir_acs,'*despiked.txt'))) #Create list of all *despiked.dat files
    
    wl = pd.read_csv('wavelength.txt', header=None, sep='\t')
    n_wl = len(wl)
    
    data_a = pd.DataFrame()
    data_c = pd.DataFrame()
    data_t = pd.Series()
    
    try:
        
        for file in fnames:
            #Read data file
            print(f'Reading {file}')      
            data = pd.read_csv(file, header=None, sep=',')
            a = data.iloc[:,1:n_wl+1]
            c = data.iloc[:,n_wl+1:]
            t = data.iloc[:,0]
            data_a = data_a.append(a)
            data_c = data_c.append(c)
            data_t = data_t.append(t)
            del data, a, c, t
        
        data = pd.concat([data_t, data_a, data_c], axis = 1)    
        data.rename(columns = {0:'Date/Time'}, inplace = True) 
        data.datetime = pd.to_datetime(data.datetime)
    
        #1-min bin and sorting data in time-ascending order.
        data_resample = data.groupby(pd.Grouper(key='Date/Time', freq='1min')).median().dropna()
        data_resample = data_resample.sort_values(by=['Date/Time'], ascending=True)
        data_resample.insert(0,'Date/Time',data_resample.index)
    
        pd.DataFrame(data=data_resample).to_csv('acs_extracted_despiked_merged_1min.txt', 
                    index=False, header=False, encoding='utf-8')
        
        print('Data merged and 1-minute binned!'+'\n')
        
    except FileNotFoundError:
        print('Error: Unable to open file: *despiked.txt!')


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------      
def rd_pangae_tempsal_underway(tsg_pangae):
    
    '''
    Thi function reads temperature and salinity data downloaded from Pangae 
    and linearly interpolates them into 1-minute interval.
    
    Input:
    tsg_pangae - absolute path of the TS data file "*oce.tab".
    
    Output - extracted and interpolated data saved as "*pyinput.txt", 
    in which the columns are datetime, temperature and salinity, respectively.
    '''
    
    foo = tsg_pangae.split('/')
    print(f'Importing {foo[-1]}')
      
    try:
        
        with open(tsg_pangae, 'r') as f:
            for num, line in enumerate(f,1):
                if '*/' in line:
                    break  
            columns = f.readlines()[0]
            
        columns = columns.split('\t')
        
        data = pd.read_csv(tsg_pangae, header=None, names=columns, 
                           skiprows=num+1, sep='\t') 
        data['Date/Time'] = pd.to_datetime(data['Date/Time'])
        data.index = data['Date/Time']
        
        data_upsampled = data.resample('1min')
        data_interpolated = data_upsampled.interpolate(method='linear')
        data_interpolated['Date/Time'] = data_interpolated.index
        data_interpolated = data_interpolated.drop(['Latitude', 'Longitude', 
                                                    'Depth water [m]'], axis=1)
            
        pd.DataFrame(data=data_interpolated).to_csv(
                tsg_pangae.replace('.tab','_pyinput.txt'), index=False, 
                header=False, encoding='utf-8')
        
        print('Temperature and salinity data extraction and interpolation Finished!')  
    
    
    except FileNotFoundError:
        print('Error: Unable to open tsg_pangae file!')

        
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------      
def acs_tempsalcorr(filename, dev, tsg, tscoeff):
    
    '''
    This function performs temperature and salinity correction for ac-s data.
    
    Input:
    filename - absolute path of ac-s file named as "acs_extracted_despiked_merged_1min.txt".
    dev - absolute path of ac-s device file (*.dev).
    tsg - absolute path of the temperature and salinity data file (*.txt). 
    Inside, first column: date time in the format of "yyyy-mm-dd HH:MM:SS"; 
    second column: temperature in degrees Celsius;
    third column: salinity. Separators between columns are commas(",").
    tscoeff - absolute path of temperature and salinity correction coefficients 
          (Sullivan et al, 2006) file (*.xls).

    Output:
    Temperature and salinity corrected ac-s data saved as 
    "acs_extracted_despiked_merged_1min_tscorr.txt", in which the columns 
    are datetime, a and c, respectively.
    '''
    
    foo = filename.split('/')
    print(f'Temperature and salinity correction for {foo[-1]} is going')
       
    #read temparature value during factory calibration (tcal) from device file.
    #(salinity value during factory calibration is 0, i.e. scal=0)
    with open(dev, 'r') as f:
        for num, line in enumerate(f,1):
            if 'tcal:' in line:
                break  
        line = line.split(' ')
        tcal = float(line[1])
    
    wl = pd.read_csv('wavelength.txt', header=None, sep='\t')   
    
    #temperature & salinity dependency coefficients (Sullivan et al, 2006).
    tscoeff = pd.read_excel(tscoeff, index_col=None, header=None, skiprows=1) 
    func_temp = inp.interp1d(tscoeff.iloc[:,0], tscoeff.iloc[:,1], 
                        'linear', fill_value='extrapolate')
    psi_temp = func_temp(wl)

    func_sala = inp.interp1d(tscoeff.iloc[:,0], tscoeff.iloc[:,5], 
                        'linear', fill_value='extrapolate')
    psi_sala = func_sala(wl)
    
    func_salc = inp.interp1d(tscoeff.iloc[:,0], tscoeff.iloc[:,3], 
                        'linear', fill_value='extrapolate')
    psi_salc = func_salc(wl)
              
    #match ac-s data with TS data
    data_acs = pd.read_csv(filename, header=None, sep=',') 
    data_acs.rename(columns = {0:'Date/Time'}, inplace = True) 
    data_tsg = pd.read_csv(tsg, header=None, sep=',') 
    data_tsg.rename(columns = {0:'Date/Time'}, inplace = True) 
        
    data_overlap = pd.merge(data_acs, data_tsg , on='Date/Time', how='inner',
                            suffixes=[" ", "_tsg"])
    a_overlap = data_overlap.iloc[:, 1:len(wl)+1]
    a_overlap.index = range(len(a_overlap))
    a_overlap.columns = range(a_overlap.shape[1])
    
    c_overlap = data_overlap.iloc[:,len(wl)+1:2*len(wl)+1]
    c_overlap.index = range(len(c_overlap))
    c_overlap.columns = range(c_overlap.shape[1])
    
    temp_delta = data_overlap.iloc[:,-2] - tcal
    
    a_tscorr = pd.DataFrame()
    c_tscorr = pd.DataFrame()
    for row in range(len(data_overlap)):
        tmp_a = a_overlap.values[row,:][:,np.newaxis] - temp_delta.iloc[
                row] * psi_temp - data_overlap.iloc[row, -1] * psi_sala
        a_tscorr = a_tscorr.append(Series(tmp_a.reshape(len(wl))), 
                                   ignore_index=True)
        
        tmp_c = c_overlap.values[row,:][:,np.newaxis] - temp_delta.iloc[
                row] * psi_temp - data_overlap.iloc[row, -1] * psi_salc
        c_tscorr = c_tscorr.append(Series(tmp_c.reshape(len(wl))), 
                                   ignore_index=True)
        del tmp_a, tmp_c
        
    
    #Write output
    output = np.hstack([data_overlap['Date/Time'].values[:,np.newaxis],
                        a_tscorr.values, c_tscorr.values])
            
    pd.DataFrame(data=output, columns=range(2*len(wl)+1)).to_csv(
            filename.replace('.txt','_tscorr.txt'), index=False, 
            header=False, encoding='utf-8')
        
    print(f'Temperature and salinity correction for {foo[-1]} Finished!')  
    

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------      
def acs_calc_apcp(filename, minute1, minute2, af_cvThold, cf_cvThold, 
                  equality=None):
    
    '''
    This function calculates ap (particulate absorption coefficient) and cp
    (particulate beam attenuation coefficient) from ac-s inline system..
    
    Input:
    filename - absolute path of ac-s file named as 
    "acs_extracted_despiked_merged_1min_tscorr.txt".
    minute1 - starting minute of stablized filtering period. e.g. 40
    minute2 - ending minute of stablized filtering period. e.g. 45. 
    minute1 < minute2.
    a(c)f_cvThold - threshold for the median coefficient of variation of 
    a(c) measurements during filtering period defined by the parameters "minute1"
    and "minute2". ac data measured during this period whose relative median 
    coefficient of variations greater than a(c)f_cvThold are considered unstable
    and removed. e.g. 0.2 or 0.5.
    equality - if equality=1, c measurements during filtering period are set 
    to be the same as a measurements (ignoring instrumental drift); 
    if equality=None (by default), a and c measurements during filtering period are "true"
    measured values and independent on each other.

    Output:
    The a and c data of particulate and dissovled materials saved as
    "acs_p.txt" and "acs_filter.txt", respectively, in which the columns 
    are datetime, a and c, respectively.
    '''
    
    foo = filename.split('/')
    print(f'Calculating ap, cp, a_filter, c_filter from {foo[-1]} ...')

    minute1 = datetime.timedelta(minutes=minute1)
    minute2 = datetime.timedelta(minutes=minute2)
    
    wl = pd.read_csv('wavelength.txt', header=None, sep='\t')
    n_wl = len(wl)
    pos_NIR = np.where(wl>=700)[0]

    data_acs = pd.read_csv(filename, header=None, sep=',')    
    t = pd.to_datetime(data_acs.iloc[:,0])
    a = data_acs.iloc[:,1:n_wl+1]
    c = data_acs.iloc[:,n_wl+1:2*n_wl+1]
    
    date_hour = pd.to_datetime(data_acs[0])
    date_hour = date_hour.apply(lambda x: x.strftime('%Y-%m-%d %H:00:00'))
    unique_hour = pd.to_datetime(date_hour.unique())
    
    af_median = pd.DataFrame()
    t_af = pd.DataFrame()
    rsd_af = pd.DataFrame()
    cf_median = pd.DataFrame()
    t_cf = pd.DataFrame()
    rsd_cf = pd.DataFrame()
   
    for i in range(len(unique_hour)):
        idx, tmp_af_rsd, tmp_af = None, None, None
        idx = np.where((t >= unique_hour[i] + minute1) & 
                       (t <= unique_hour[i] + minute2))[0]
        if len(idx)>0:
            tmp_af_rsd = rsd_median(a.iloc[idx,:pos_NIR[0]-1])

            tmp_cf_rsd= rsd_median(c.iloc[idx,:pos_NIR[0]-1])
            
            if any(tmp_af_rsd > af_cvThold):
                a.iloc[idx,:] = np.nan
            if any(tmp_cf_rsd > cf_cvThold):
                c.iloc[idx,:] = np.nan
            
            tmp_af = np.nanmedian(a.iloc[idx,:], axis=0)
            tmp_cf = np.nanmedian(c.iloc[idx,:], axis=0)
            
            if all(np.isfinite(tmp_af)):
                tmp_ta = pd.Series(unique_hour[i] + minute1)
                t_af = t_af.append(tmp_ta, ignore_index=True)
                tmp_af = tmp_af[np.newaxis,:]
                af_median = af_median.append(pd.DataFrame(tmp_af), 
                                             ignore_index=True)
                tmp_af_rsd = tmp_af_rsd[np.newaxis,:]
                rsd_af = rsd_af.append(pd.DataFrame(tmp_af_rsd), 
                                       ignore_index=True)
                
            if all(np.isfinite(tmp_cf)):
                tmp_tc = pd.Series(unique_hour[i] + minute1)
                t_cf = t_cf.append(tmp_tc, ignore_index=True)
                tmp_cf = tmp_cf[np.newaxis,:]
                cf_median = cf_median.append(pd.DataFrame(tmp_cf), 
                                             ignore_index=True)
                tmp_cf_rsd = tmp_cf_rsd[np.newaxis,:]
                rsd_cf = rsd_cf.append(pd.DataFrame(tmp_cf_rsd), 
                                       ignore_index=True)    

    af_median.insert(0, 'Date/Time', t_af.squeeze())
    t= t.to_frame('Date/Time')
    af_interp = pd.merge(t,af_median, on='Date/Time', 
                                   how='outer').sort_values('Date/Time')
    af_interp.index = af_interp['Date/Time']
    af_interp = af_interp.drop('Date/Time', axis=1)
    af_interp = af_interp.interpolate(method='linear', 
              limit_direction ='both')
    
    if equality==1:
        cf_interp = af_interp
    else:
        cf_median.insert(0, 'Date/Time', t_cf.squeeze())
        cf_interp = pd.merge(t,cf_median, on='Date/Time', 
                                   how='outer').sort_values('Date/Time')
        cf_interp.index = cf_interp['Date/Time']
        cf_interp = cf_interp.drop('Date/Time', axis=1)
        cf_interp = cf_interp.interpolate(method='linear', 
                                          limit_direction ='both')
    
    ap = a.values - af_interp.values
    cp = c.values - cf_interp.values
     
    t_minutes = t.squeeze().apply(lambda x: datetime.timedelta(minutes=x.minute))
    minutes_filter = np.where((t_minutes >= minute1) & (t_minutes <= minute2))[0]
    ap = np.delete(ap, minutes_filter, axis = 0)
    cp = np.delete(cp, minutes_filter, axis = 0)
    tp = np.delete(data_acs.iloc[:,0].values[:, np.newaxis], minutes_filter,
                   axis = 0)
       
    #Write ac-s output during filtering periods
    output_filter = np.hstack([data_acs.iloc[:,0].values[:, np.newaxis], 
                               af_interp.values, cf_interp.values])
            
    pd.DataFrame(data=output_filter, columns=range(2*n_wl+1)).to_csv(
            'acs_filter.txt', index=False, header=False, encoding='utf-8')
    
    #Write ac-s output (i.e. ap & cp) during non-filtering periods
    output_p = np.hstack([tp, ap, cp])
    pd.DataFrame(data=output_p, columns=range(2*n_wl+1)).to_csv(
            'acs_p.txt', index=False, header=False, encoding='utf-8')

    print(f'Calculation of ap, cp, a_filter, c_filter from {foo[-1]} Finished!')  


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def costf_tscat_slade(delta_temp, ap, cp, psi_temp, pos_NIR, pos_ref, 
                      scatcorr=None):
    
    '''
    This function defines the cost function according to or adapted from
    Equation 6 in Slade et al. (2010) for residual temperature correction for 
    cp data and combined residual temperature and scattering correction for 
    ap data. 
    
    Input:
    Data type of all input arguments and parameters are numpy array.
    delta_temp - residual temperature value between water samples and TSG 
    measured temperature.
    a(c)p - TS corrected (Sullivan et al. 2006) particulate absorption 
    (attenuation) coefficients.
    psi_temp - temperature correction coefficients from Sullivan et al. 2006.
    pos_NIR - indexes of near infrared wavelengths for ap/cp spectra.
    pos_ref - index of reference wavelength for scattering correction of 
    absorption measurements (ac-9: 715 nm, ac-s: 730 nm).
    scatcorr - Scattering correction method. If scatcorr=None (by default), 
    it takes the proportional method from Zaneveld et al (1994); 
    if scatcorr=1, scattering correction follows Roettgers et al (2013).
    
    Output:
    Cost function according to Equation 6 in Slade et al. (2010).
    '''
    
    bp = cp - ap
    
    if scatcorr==1: 
        #Zaneveld et al. (1994)
        costf = sum( abs( ap[pos_NIR] - 
                     psi_temp[pos_NIR] * delta_temp
                     - (ap[pos_ref]
                     - psi_temp[pos_ref] * delta_temp) ) )
    else: 
        #Roettgers et al. (2013)
        costf = sum( abs( ap[pos_NIR] - 
                     psi_temp[pos_NIR] * delta_temp
                     - ( (ap[pos_ref]
                     - psi_temp[pos_ref] * delta_temp)
                     /bp[pos_ref] ) * bp[pos_NIR] ) )
        
    return costf


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------  
def acs_residtemp_scatcorr(filename, tscoeff, scatcorr=None):
    
    '''
    This function corrects the residual temperature and scattering errors of 
    ap and residual temperature of cp from ac-s inline system. 
    Residual temperature correction of ap and cp follows Slade et al. (2010).
    
    Input:
    filename - absolute path of ac-s file named as "acs_p.txt".
    tscoeff - absolute path of temperature and salinity correction coefficients 
          (Sullivan et al, 2006) file (*.xls).
    scatcorr - Scattering correction method. If scatcorr=None (by default), 
    it takes the proportional method from Zaneveld et al. (1994); 
    if scatcorr=1, scattering correction follows Roettgers et al. (2013).
    
    Output:
    The finally quality controlled a and c data of particulate materials saved as
    "acs_TSalScatCorr.txt", in which the columns are datetime, a and c, respectively.
    '''
    
    foo = filename.split('/')
    print(f'residual temperature and scattering correction for \
          {foo[-1]} is going...')
    
    #wavelengths
    wl = pd.read_csv('wavelength.txt', header=None, sep='\t')
    n_wl = len(wl)
    pos_NIR = np.where(wl>=700)[0]
    pos_ref = np.where(wl>=730)[0][0]
    pos_440 = np.where(wl>=440)[0][0]
    pos_675 = np.where(wl>=675)[0][0]
    
    #temperature correction coefficients (Sullivan et al. 2006).
    tscoeff = pd.read_excel(tscoeff, index_col=None, header=None, skiprows=1) 
    func_temp = inp.interp1d(tscoeff.iloc[:,0], tscoeff.iloc[:,1], 
                        'linear', fill_value='extrapolate')
    psi_temp = func_temp(wl).squeeze()

    #ac-s data
    data_acs = pd.read_csv(filename, header=None, sep=',')    
    ap = data_acs.to_numpy()[:,1:n_wl+1].astype('float')
    cp = data_acs.to_numpy()[:,n_wl+1:2*n_wl+1].astype('float')
    tp = data_acs.iloc[:,0].reindex(range(len(ap)))
    
    
    ap_TSalScatCorr = np.empty((len(ap), n_wl)) * np.nan
    cp_TSalCorr = np.empty((len(ap), n_wl)) * np.nan
    delta_temp = np.zeros(len(ap))* np.nan
    fiterr = np.zeros(len(ap))* np.nan
    
    
    for i in range(len(ap)):
        if all(np.isfinite(ap[i,:])):
            delta_temp0 = 0
            # minimization routine (Equation 6 in Slade et al. 2010)
            delta_temp[i] = scipy.optimize.fmin(func=costf_tscat_slade, 
                                x0=delta_temp0, 
                                args=(ap[i,:], 
                                      cp[i,:], 
                                      psi_temp, 
                                      pos_NIR, 
                                      pos_ref,
                                      scatcorr), 
                                      xtol=1e-8, ftol=1e-8, 
                                      maxiter=20000000, maxfun=20000,
                                      disp=0)
 
            fiterr[i] = costf_tscat_slade(delta_temp[i], ap[i,:], 
                                      cp[i,:], 
                                      psi_temp, 
                                      pos_NIR, 
                                      pos_ref,
                                      scatcorr)

            bp = cp[i,:] - ap[i,:]
            
            if scatcorr==1:
                #scattering correction following Roettgers et al. (2013).
                ap_residtemp_ref = ap[i,pos_ref] - psi_temp[pos_ref]*delta_temp[i]
                
                ap_ref = ap_residtemp_ref - 0.212 * math.copysign(
                        1, ap_residtemp_ref) * abs(ap_residtemp_ref)**1.135
                
                ap_TSalScatCorr[i,:] = ap[i,:] - psi_temp*delta_temp[i] - (
                        ap_ref/bp[pos_ref]) * bp
                   
            else:
                #scattering correction following Zaneveld et al. (1994).
                ap_TSalScatCorr[i,:] = ap[i,:] - psi_temp*delta_temp[i] - (
                    (ap[i,pos_ref] - psi_temp[pos_ref]*delta_temp[i])
                    /bp[pos_ref])*bp
            
            
            cp_TSalCorr[i,:] = cp[i,:] - psi_temp*delta_temp[i]
        
        
    bp = cp_TSalCorr - ap_TSalScatCorr
    
    pos1 = [ i for i in range(len(ap_TSalScatCorr)) \
            if any(ap_TSalScatCorr[i,1:pos_NIR[0]-1]<=0) ]
    pos2 = [ i for i in range(len(ap_TSalScatCorr)) \
            if ap_TSalScatCorr[i,pos_440]/ap_TSalScatCorr[i,pos_675]<=0.95 ]
    pos3 = [ i for i in range(len(ap_TSalScatCorr)) \
            if any(bp[i,1:pos_NIR[0]-1]<=0) ]
    pos4 = [ i for i in range(len(ap_TSalScatCorr)) \
            if any(ap_TSalScatCorr[i,1:5]<=0.001) ]
    pos5 = [ i for i in range(len(ap_TSalScatCorr)) \
            if np.isnan(ap_TSalScatCorr[i,1])]
    pos6 = [ i for i in range(len(ap_TSalScatCorr)) \
            if np.isnan(cp_TSalCorr[i,1])]

    pos = list(set(pos1) |set(pos2) |set(pos3) |set(pos4) |set(pos5) |set(pos6))
    ap_TSalScatCorr = np.delete(ap_TSalScatCorr, pos, axis = 0)
    cp_TSalCorr = np.delete(cp_TSalCorr, pos, axis = 0)
    bp = np.delete(bp, pos, axis = 0)
    delta_temp = np.delete(delta_temp, pos, axis = 0)
    fiterr = np.delete(fiterr, pos, axis = 0)
    tp = tp.drop(index=pos)
    
    
    #Write ac-s output 
    output= np.hstack([tp.values[:,np.newaxis], ap_TSalScatCorr.astype('object'),\
                       cp_TSalCorr.astype('object')])
            
    pd.DataFrame(data=output, columns=range(2*n_wl+1)).to_csv(
            'acs_TSalScatCorr.txt', index=False, header=False, encoding='utf-8')

            
    print(f'Residual temperature correction for cp data and combined \
          residual temperature and scattering correction for ap data for \
          {foo[-1]} Finished!')
    

#-----------------------------------------------------------------------------
# Main Program. Call functions here.
#-----------------------------------------------------------------------------
if __name__ == '__main__':
    
    #directory of ac-s raw data. Modify path!
    dir_acs = '/Users/yliu/Data/acs_underway_py/PS113_acs213/' 
    #path of temperature and salinity data file downloaded from pangae. optional. Modify path!
    tsg_pangae='/Users/yliu/Data/acs_underway_py/TSG/PS113_surf_oce.tab'
    #path of extracted temperature and salinity file. Modify path!
    tsg='/Users/yliu/Data/acs_underway_py/TSG/PS113_surf_oce_pyinput.txt'
    #path of device file. Modify path!
    dev='/Users/yliu/Data/acs_underway_py/device/acs213_092716.dev'
    #path of the temperature and salinity correction coefficients file. Modify path!
    tscoeff='/Users/yliu/Data/acs_underway_py/TSG/Sullivan_etal_2006_instrumentspecific.xls'

    #Set working directory
    os.chdir(dir_acs) 
    print ('Switching to working directory:'+ str(os.getcwd()))
    
    #Extract raw ac-s data from wetview or compass output file by file.    
    fnames_raw=sorted(glob.glob(os.path.join(dir_acs,'*.dat')))
    for file in fnames_raw:
#        acs_rd_wetview(file)
        acs_rd_compass(file)
    
    #De-spike raw ac-s data file by file. Modify function arguments!
    fnames_extracted=sorted(glob.glob(os.path.join(dir_acs,'*extracted.txt')))
    for file in fnames_extracted:
        acs_rm_spikes(file, a_Thold=1, c_Thold=2, a_cvThold=5, c_cvThold=10)
    
    #Merge all raw ac-s files and bin data into 1-minute interval.
    acs_filemerge_oneminbin(dir_acs)
    
    #Extract temperature and salinity data from TS file downloaded from pangae.
    rd_pangae_tempsal_underway(tsg_pangae)
    
    #Remperature and salinity correction. Modify path!
    filename = '/Users/yliu/Data/acs_underway_py/PS113_acs213/acs_extracted_despiked_merged_1min.txt'
    acs_tempsalcorr(filename, dev, tsg, tscoeff)
    
    #Calculate particulate absorption and attenuation. Modify path and function arguments!
    filename = '/Users/yliu/Data/acs_underway_py/PS113_acs213/acs_extracted_despiked_merged_1min_tscorr.txt'
    acs_calc_apcp(filename, minute1=41, miute2=45, af_cvThold=0.5, cf_cvThold=0.5, equality=1)
    
    #Residual temperature and scattering correction. Modify path and function arguments!
    filename = '/Users/yliu/Data/acs_underway_py/PS113_acs213/acs_p.txt'
    acs_residtemp_scatcorr(filename, tscoeff, scatcorr=1)

