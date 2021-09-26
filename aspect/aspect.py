import os
import arcpy
import pandas

df= pandas.read_csv('C:/Users/Daniel/Documents/Research/O2_basins.csv')
folder = "C:/Users/Daniel/Documents/Research/gis_projects/O2_catchments/O2_rasters"

os.chdir(folder)
arcpy.env.workspace = folder

arcpy.CheckOutExtension("spatial")

#catchment = "07208500" #USGS gauge code
#dem_path = "E_" + catchment

sitelist=df['SITE_NO_tx'].tolist()
for i in range(0,len(sitelist)):
    sitelist[i]=str(sitelist[i])
    if len(sitelist[i])<8:
        sitelist[i]= '0'+sitelist[i]

x = 0
sitenum = sitelist[x]
Eraster = folder+"/E_"+sitenum+'.tif'
Araster = folder+"/A_"+sitenum+'.tif'
Hraster = folder+"/H_"+sitenum+'.tif'


outAspect = arcpy.ddd.Aspect(Eraster, Araster)

outHS = arcpy.ddd.HillShade(Eraster, Hraster, 180, 37, "NO_SHADOWS", 1)

