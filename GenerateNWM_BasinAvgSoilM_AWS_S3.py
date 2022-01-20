# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 08:26:19 2022

The code reads noaa-nwm volumetric soil moisture from AWS S3 buckets (without saving data)
and generate basin average soilmoisture for each basin

input:
    1. noaa-aws aws s3
    2. weight tables
    3. basins shapefile
    
output:
    1. seperate csv file for each basins including soil moisture basin avg, 
    each csv file includes timeseries for all the four layers
    
    
*** NOTE: it will be good to add a loop to run the code monthly or yearly and 
        Save the data for each month/year seperately***

@author: sardekani
"""

import boto3
import s3fs
import xarray as xr
import zarr
import numpy as np
import pandas as pd
import geopandas as gpd
from datetime import datetime, timedelta
import timeit
import os

print('NOTE: it will be good to add a loop to run the code monthly or yearly and Save the data for each month/year seperately')

inputDir="R:/2022/LCRA/NWM_SM/BasinMasks/BasinWeightTable/"
outdir='R:/2022/LCRA/NWM_SM/BasinAvg3/'
if not os.path.exists(outdir): os.mkdir(outdir)

# load the shapefile includings all basins
Shapefile="P:/2020/LCRA/BasinDelineationLKC_Jan8/SubWatershedSHP/SWAllMergedDis_NWM.shp"
gdb=gpd.read_file(Shapefile)
id_no_Ar=list(gdb.ProjID)

# read noaa-nws volumetric soil moisture data from s3 aws
path='https://noaa-nwm-retrospective-2-1-zarr-pds.s3.amazonaws.com/ldasout.zarr'
ds = xr.open_zarr(path,consolidated=True)
ds_SOIL=ds.SOIL_W

# clean up the shop a little bit
del(ds)

# Start date and end date of data
d1=datetime.utcfromtimestamp(ds_SOIL.time[0].values.tolist()/1e9)
d2=datetime.utcfromtimestamp(ds_SOIL.time[-1].values.tolist()/1e9)
d2=datetime(1979, 3, 1, 3, 0)
# Generate date timeseries
dateAr=pd.date_range(start=d1, end=d2, freq='3H')

d1_str=datetime.strftime(d1, '%Y%m%d%H')
d2_str=datetime.strftime(d1, '%Y%m%d%H')

# Find the minimum and maximum col and row values
cminAr=[]
cmaxAr=[]
rminAr=[]
rmaxAr=[]
for index, id_no in enumerate(id_no_Ar):
    infile=inputDir+str(id_no)[2:]+'_NWM_Weight.csv'
    df=pd.read_csv(infile)
    cminAr.append(df.col.min())
    cmaxAr.append(df.col.max())
    rminAr.append(df.row.min())
    rmaxAr.append(df.row.max())
cmin=min(cminAr)
cmax=max(cmaxAr)+1
rmin=min(rminAr)
rmax=max(rmaxAr)+1

print(cmin,cmax,rmin,rmax)

start = timeit.default_timer()

# clean up the VSM_crop to make sure it is not in the memory for later
try: del(VSM_crop)
except: print('VSM_crop is not in the memory')

# find the index of start and end date time
d_ind1=int((dateAr[0]-d1)/ np.timedelta64(3, 'h'))
d_ind2=int((dateAr[-1]-d1)/ np.timedelta64(3, 'h'))+1

# iterate through each basin and calculate basin average timeseries
for index, id_no in enumerate(id_no_Ar):
    print(index)
    # read weigh table 
    infile=inputDir+str(id_no)[2:]+'_NWM_Weight.csv'
    df=pd.read_csv(infile)
    col_l=df.col
    row_l=df.row
    w=df.weight
    # create empty lists to write basin average timeseries in later
    SOIL_W1=[]
    SOIL_W2=[]
    SOIL_W3=[]
    SOIL_W4=[]
    
    # Crop the data to reduce the process time
    try: VSM_crop
    except:  VSM_crop=ds_SOIL[d_ind1:d_ind2,cmin:cmax,:,rmin:rmax].values

    for i in range(len(dateAr)):
        # create empty lists to later write volumetric soil moisture values in
        VSM_Layer0=[]
        VSM_Layer1=[]
        VSM_Layer2=[]
        VSM_Layer3=[]
        for j in range(len(df)):
            # extract SW 
            col=col_l[j]-cmin
            row=row_l[j]-rmin
            VSM=VSM_crop[i,col,:,row].tolist()
            VSM_Layer0.append(VSM[0])
            VSM_Layer1.append(VSM[1])
            VSM_Layer2.append(VSM[2])
            VSM_Layer3.append(VSM[3])

           
        # calculate Soil moisture basin avg 
        SOIL_W1.append(np.nansum(VSM_Layer0*w)/sum(w))
        SOIL_W2.append(np.nansum(VSM_Layer1*w)/sum(w))
        SOIL_W3.append(np.nansum(VSM_Layer2*w)/sum(w))
        SOIL_W4.append(np.nansum(VSM_Layer3*w)/sum(w))
        
        
    # Write data into a csv file for each basin
    df_SOIL=pd.DataFrame(zip(dateAr.tolist(),SOIL_W1,SOIL_W2,SOIL_W3,SOIL_W4))
    df_SOIL.columns=['date','SOIL_W1','SOIL_W2','SOIL_W3','SOIL_W4']
    outfile=outdir+d1_str+'_'+str(id_no)[2:]+'_BasinAvg.csv'
    df_SOIL.to_csv(outfile, index=False)
    
stop = timeit.default_timer()
print('Time: ', stop - start) 