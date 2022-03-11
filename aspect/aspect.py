

print("importing libraries...")

import os
import arcpy
import pandas

print("reading O2_basins.csv...")
df= pandas.read_csv('C:/Users/Daniel/Documents/Research/O2_basins.csv')
folder = "C:/Users/Daniel/Documents/Research/gis_projects/O2_catchments/O2_rasters"

os.chdir(folder)
arcpy.env.workspace = folder

arcpy.CheckOutExtension("spatial")

#catchment = "07208500" #USGS gauge code
#dem_path = "E_" + catchment

print('formatting site codes...')
sitelist=df['SITE_NO_tx'].tolist()
for i in range(0,len(sitelist)):
    sitelist[i]=str(sitelist[i])
    if len(sitelist[i])<8:
        sitelist[i]= '0'+sitelist[i]

x = input("enter feature ID for desired catchment:") #User input: Enter the FID from ArcGIS of the catchment you want to process.
sitenum = sitelist[x]
print "processing DEM for ", df['STANAME'][x]
Eraster = folder+"/E_"+sitenum+'.tif'
Araster = folder+"/A_"+sitenum+'.tif'
Hraster = folder+"/H_"+sitenum+'.tif'

print('generating aspect raster...')
outAspect = arcpy.ddd.Aspect(Eraster, Araster)

print('generating hillshade raster...')
outHS = arcpy.ddd.HillShade(Eraster, Hraster, 180, 37, "NO_SHADOWS", 1)

print('generating slope raster...')

print('all rasters complete')
