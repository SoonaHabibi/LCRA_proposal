# -*- coding: utf-8 -*-
"""
Created on Tue Oct 06 11:05:54 2015

@author: lcunha
"""

# Import arcpy module
import csv

import sys 
sys.path.append("C:/Program Files (x86)/ArcGIS/Desktop10.1/arcpy");
sys.path.append("C:/Python27/ArcGIS10.1/Lib");
sys.path.append("C:/Program Files (x86)/ArcGIS/Desktop10.1/bin");
sys.path.append("O:/Python/WEST_Python_FunctionAug2019/")
import subprocess
import BasicFunction_py3 as BF   
import numpy as np
import pandas as pd
#import arcpy
import numpy as np
import string 
import os
from pylab import *
import matplotlib
import math as m
from math import cos, sin, asin, sqrt, radians
#from arcpy import env
#from arcpy.sa import *
import glob
#import arcgisscripting
import os

from osgeo import gdal, ogr
from osgeo import osr

print ("TODO: Does not WORK")
OutputDir="R:/2022/LCRA/NWM_SM/BasinMasks/BasinWeightTable/"
if not os.path.exists(OutputDir):    os.mkdir(OutputDir) 
OutputDirSHP="R:/2022/LCRA/NWM_SM/BasinMasks/SHP/"
if not os.path.exists(OutputDirSHP):    os.mkdir(OutputDirSHP) 

# Save the tif mask file if Flag==1 
Flag=0 

Name="_NLDAS";Filename = 'O:/NLDAS/NLDAS_Output/NOAH/NLDAS_SOILM_19790102.0600SOILM_0-80mm.tif';ntimes=26
Name="_AORC";Filename = 'O:/Rainfall_MPE_AORC/WGRFC/WGRFC_Mask.tif';ntimes=10
Name="_NWM";Filename = '//westfolsom/D/WESTDataCode/NWS_NationalWaterModel/CompFiles_v2.0/Tiff/SOIL_W199912300600.tif';ntimes=10 ;projtif='NWM'

img = gdal.Open(Filename)
band = img.GetRasterBand(1)
ncol=img.RasterXSize 
nrow=img.RasterYSize
geotransform = img.GetGeoTransform()
leftx=geotransform[0]
uppery=geotransform[3]
resdegr=geotransform[1]
ProjTif=img.GetProjection()


MaskData=pd.read_csv("P:/2020/LCRA/BasinDelineationLKC_Jan8/SummaryResWatAllSitesJustUpTTCalculation_2498fix.csv") 

#MaskData=pd.read_csv("P:/2019/NASA_UTA/USACE/AORC_Rain/SummaryResWatAllSitesJustUp.csv") 
Shapefile="P:/2020/LCRA/BasinDelineationLKC_Jan8/SubWatershedSHP/SWAllMergedDis_"+projtif+".shp"
source_ds = ogr.Open(Shapefile)
source_layer = source_ds.GetLayer()
#GeogCS = geom.GetSpatialReference()
spatialRef = source_layer.GetSpatialRef()
spatialRefWkt = spatialRef.ExportToWkt()
nr_features = source_layer.GetFeatureCount()
proj = osr.SpatialReference(wkt=spatialRefWkt)

id_no_Ar=[]
# Generate separate shapefiles for each feature of the main shapefile
for feat_id in range(len(MaskData)):
    feature = source_layer.GetFeature(feat_id)
    transformed = feature.GetGeometryRef()
    id_no = feature.GetField("ProjID")                        # Field name could be different based on the shapefile
    print (str(feat_id) + " " +str(id_no) )
    id_no_Ar.append(id_no)
    
    outShapefile = os.path.join(OutputDirSHP, id_no+"_"+projtif+".shp")
    outDriver = ogr.GetDriverByName("ESRI Shapefile")    
    
    # Remove output shapefile if it already exists
    if os.path.exists(outShapefile): outDriver.DeleteDataSource(outShapefile)        
    
    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(outShapefile)                        
    outLayer = outDataSource.CreateLayer(id_no, geom_type=ogr.wkbPolygon,srs=spatialRef)
    
    # Add a presence field
    presField = ogr.FieldDefn("ProjID", ogr.OFTString)                                
    outLayer.CreateField(presField)                            # Field name could be different based on the shapefile

    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    feat = ogr.Feature(featureDefn)
    geom = feature.GetGeometryRef()

    feat.SetGeometry(geom)
    feat.SetField("ProjID", id_no)
    outLayer.CreateFeature(feat)
    feat = None

    # Save and close DataSource
    outDataSource = None
 
    
for i in range(len(MaskData)):
    
    id_no = MaskData['LCRAID'].iloc[i]
    Shapefile = OutputDirSHP+'SW'+ str(id_no)+"_"+projtif+".shp"

    print (Shapefile)    


    if (not os.path.exists(Shapefile)): 
        print ("Shapefile does not exist:" +Shapefile)
    else:
        print ("Shapefile exist: " +Shapefile)    
               
        RasterFile=OutputDir+str(id_no)+"MP_HR"+"_"+projtif+".tif"
        if os.path.exists(RasterFile): os.remove(RasterFile)
        if (not os.path.exists(RasterFile)): 
            print( 'Not Available: ', RasterFile)
            geotransformH=np.copy(geotransform)
            geotransformH[1]=geotransform[1]/ntimes
            geotransformH[5]=-geotransform[1]/ntimes
            nrowH=nrow*ntimes
            ncolH=ncol*ntimes
            resdegrH=resdegr/ntimes
            
            #BF.CreateMaskFromShapefile(Shapefile,RasterFile,geotransformH,spatialRef,nrowH,ncolH,resdegrH, leftx, uppery)
            source_ds2 = ogr.Open(Shapefile)
            source_layer2 = source_ds2.GetLayer()
        
            driver = gdal.GetDriverByName('GTiff')
            
            #target_ds = driver.Create(RasterFile, int(ncol), int(nrow), 1, gdal.GDT_Float32 )
            target_ds = driver.Create(RasterFile, int(ncolH), int(nrowH), 1, gdal.GDT_Float32 )    
            target_ds.SetGeoTransform(geotransformH)       
            srs = osr.SpatialReference()
            #srs.SetWellKnownGeogCS(GeogCS )
            target_ds.SetProjection(spatialRef.ExportToWkt() )
            res=target_ds.SetProjection( ProjTif )
            if res != 0:
                print('Setting projection failed ' + str(res))            
            
            #gdal.RasterizeLayer(target_ds, [1], source_layer, None, None, burn_values=[1], ['ALL_TOUCHED=TRUE'])
            gdal.RasterizeLayer(target_ds, [1], source_layer2, None, None, [1], ['ALL_TOUCHED=TRUE'])
            target_ds = None

    
        img = gdal.Open(RasterFile)
        ncolH=img.RasterXSize 
        nrowH=img.RasterYSize
        geotransformH = img.GetGeoTransform()    
        resdegrH=geotransformH[1]
        band = img.GetRasterBand(1)
        dataH = band.ReadAsArray()
           
        dataRef=np.zeros((int(nrowH),int(ncolH)))
        dataRef[dataH>0]=1
        #matrixBasin=BF.CreateMaskMatrix(dataRef,geotransformH,nrow,ncol,resdegr, leftx, uppery)
        matrixBasin=BF.CreateMaskMatrix(dataRef,geotransformH,nrow,ncol,resdegr, leftx, uppery)
        
        matrixBasin=100.0*matrixBasin/float(pow(resdegr/resdegrH,2))
        matrixBasin[matrixBasin>100]=100.


        # Create weigth csv file
        size=matrixBasin.shape
        listWeight=matrixBasin.flatten().tolist()              # Convert weight matrix into a list
        row=np.repeat(list(range(size[0])),size[1])             # row's index of 'listWeight' elements in matrix basin 
        col=list(range(size[1]))*size[0]                        # column's index of 'listWeight' elements in matrix basin 
        weightDF=zip(row,col,listWeight)
        frameEvents = pd.DataFrame(weightDF,columns=['row','col','weight']) 
        w=frameEvents[frameEvents['weight']!=0]                  # remove all rows with weight = 0
        
        CSVout=OutputDir+str(id_no)+Name+'_Weight.csv'
        w.to_csv(CSVout,index=False)
        FileName=OutputDir+str(id_no)+Name+'.tif'
        print(FileName)
        
        # write mask array into tif
        if Flag==1:
            BF.CreateMatrixFileFloat(FileName,matrixBasin,ncol, nrow,geotransform,spatialRefWkt)
        
        img = None