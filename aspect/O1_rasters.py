

print("importing libraries...")

import os
import arcpy
import pandas

print("reading catchment_comparison.csv...")
df = pandas.read_csv('C:/Users/Daniel/Documents/Research/data/dataframes/LRO_site_codes.csv')
folder = "C:/Users/Daniel/Documents/Research/hypsometric_curves"

os.chdir(folder)
arcpy.env.workspace = folder

arcpy.CheckOutExtension("spatial")

#catchment = "07208500" #USGS gauge code
#dem_path = "E_" + catchment



x = input("enter feature ID for desired catchment:") #User input: Enter the 2 letter id.
site=df['site_abbr'][x-1]
print "processing DEM for ", df['site_name'][x-1]

Eraster = folder+"/"+site+'_DEM.tif'
Araster = folder+"/A_"+site+'.tif'
Hraster = folder+"/H_"+site+'.tif'
Sraster = folder+"/S_"+site+'.tif'

print('generating aspect raster...')
outAspect = arcpy.ddd.Aspect(Eraster, Araster)

print('generating hillshade raster...')
outHS = arcpy.ddd.HillShade(Eraster, Hraster, 180, 37, "NO_SHADOWS", 1)

print('generating slope raster...')
out_raster = arcpy.sa.Slope(Eraster, "DEGREE", 1, "PLANAR", "METER"); out_raster.save(Sraster)


print('all rasters complete')


arcpy.ddd.Slope("C:/Users/Daniel/Documents/Research/hypsometric_curves/SC_DEM.tif", r"C:\Users\Daniel\Documents\Research\gis_projects\Logan_Watershed\sp_cr_slope.tif", "DEGREE", 1, "PLANAR", "METER")