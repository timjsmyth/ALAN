#!/usr/bin/env python
import os
import sys
from timeit import default_timer as timer
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from itertools import product
from random import randrange, uniform
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import RegularGridInterpolator
from PIL import Image
from osgeo import gdal
from numpy import zeros, newaxis
from scipy.ndimage import zoom
from matplotlib.colors import LogNorm
from mpl_toolkits.basemap import Basemap
import argparse
import warnings
warnings.simplefilter('ignore', np.RankWarning)

def main():

   parser = argparse.ArgumentParser(description='Produce ALAN maps')

   parser.add_argument("--Zc",  action='store_true',  default=False, help="Output critical depth where ALAN drops below threshold.")
   parser.add_argument("--seabed",  action='store_true',  default=False, help="Output irradiance at the seafloor.")
   parser.add_argument("--cloudy",  action='store_true',  default=False, help="Use scaling factors for cloudy skies")
   parser.add_argument("--roi", type=int, nargs=4, help="ROI - Region of interest")
   parser.add_argument("--odir", type=str, default='', help="Output directory for netcdf imagery")
   args = parser.parse_args()

   filestr_box = ""
   cloudy = "clear"

   # prevent both the seabed and critical depth options being selected
   if args.seabed and args.Zc:
      print("Select one of either --Zc or --seabed")
      sys.exit()

   if args.cloudy:
      cloudy = "cloudy"

   # Region of interest 
   if args.roi:
      print("Region of interest selected")
      print(args.roi)
      roi_box = args.roi
      print("ROI: ", roi_box)
      min_lat = roi_box[0]
      max_lat = roi_box[1]
      min_lon = roi_box[2]
      max_lon = roi_box[3]
      filestr_box = '_'+str(min_lat)+'S_'+str(max_lat)+'N_'+str(min_lon)+'W_'+str(max_lon)+'E'
      print(filestr_box)

   if args.seabed:
      print("Calculating the ALAN irradiance at the seabed")
      if args.odir:
         outdir = args.odir
      else: 
         outdir = '/home/scratch/data/workspace/ALAN/seafloor/'
      filedesc = "-Iseafloor-"      
   
      print(" Reading in bathymetry data")
      if args.roi:
         bathyinfile=nc.Dataset('/home/scratch/data/workspace/ALAN/ancillary/bathymetry/gebco_global_0.5km_g7.nc', 'r')  
      else:
         bathyinfile=nc.Dataset('/home/scratch/data/workspace/ALAN/ancillary/bathymetry/gebco_global_4km_g7.nc', 'r')
      
      z_bathy=bathyinfile.variables['elevation'][:]
      bathy_lats=bathyinfile.variables['lat'][:]
      bathy_lons=bathyinfile.variables['lon'][:]
   
      z_bathy=np.flipud(z_bathy)
      bathy_lats=np.flipud(bathy_lats)
      print("  Re-Shape of bathymetry matrix")
      # Add a leading dimension so equivalent to the Kd files and 
      # multiply elevation by -1.0 to get positive as depth below surface
      z_bathy=-1.0*z_bathy[newaxis,...]
      print(z_bathy.shape)
      
   # Critical depth
   if args.Zc:
      print("Calculating the critical depth below which ALAN drops below threshold")
      if args.odir:
         outdir = args.odir 
      else:
         outdir = '/home/scratch/data/workspace/ALAN/Zc/'
      filedesc = "-Zcritical-"      

   # slope and intercept values for mCd/m2 -> uW/m2 
   if args.cloudy:
      # "CLOUDY" sky
      m_blue = 17.97854676
      m_green = 35.80029137
      m_red = 29.40875316
      c_blue = 57.22738163
      c_green = 18.58249682
      c_red = 23.73604658
   else:
      # "CLEAR" sky
      m_blue = 4.530851872
      m_green = 7.272195079
      m_red = 6.372282777
      c_blue = 59.57581038
      c_green = 25.091669
      c_red = 26.19770495

   # From Batnes et al. (2015) - table 5
   thresh_irr_total_uW_m2 = 0.102 # in uW/m2 converted from uE/m2/s

   # FALCHI map runs between 85N -> 60S (145 degrees)
   print("Reading in Falchi map")
   Image.MAX_IMAGE_PIXELS = 751939200
   im = gdal.Open("/home/scratch/data/workspace/ALAN/World_Atlas_2015.tif")
   falchi_map = im.ReadAsArray()
   width = im.RasterXSize
   height = im.RasterYSize
   gt = im.GetGeoTransform()
   minx = gt[0]
   miny = gt[3] + width*gt[4] + height*gt[5] 
   maxx = gt[0] + width*gt[1] + height*gt[2]
   maxy = gt[3]
   
   falchi_lats=np.zeros(int(height))*0.
   falchi_lons=np.zeros(int(width))*0.
   
   falchi_lats_orig = falchi_lats
   falchi_lons_orig = falchi_lons
   
   for ii in range(width):
      falchi_lons[ii] = ii*gt[1]+gt[0]
   
   for jj in range(height):
      falchi_lats[jj] = gt[3]+jj*gt[5]
   
   print("Shape of Falchi map")
   print(falchi_map.shape)

   datadir = '/home/scratch/data/workspace/ALAN/Kd/'
   for filename in sorted(os.listdir(datadir)):
      if filename.endswith("fv4.0.nc"): 
         # determine the filename for the netCDF and read in the file
         #print(os.path.join(datadir, filename))
         print("=======================================")
         print(" Input filename: ", filename)
         fname = os.path.join(datadir, filename)
         getmonth = filename.split('-')
         ofname = getmonth[0]+"-"+getmonth[1]+"-"+getmonth[2]+"-"+getmonth[3]+"-"+getmonth[4]+filedesc+cloudy+"-"+getmonth[6]+"-fv4.0"+filestr_box+".nc"
         print("  Output filename: ", ofname)

         infile=nc.Dataset(fname, 'r')
         kd_blue=infile.variables['kd_blue'][:]
         kd_green=infile.variables['kd_green'][:]
         kd_red=infile.variables['kd_red'][:]

         nc_dims = [dim for dim in infile.dimensions]  # list of nc dimensions
         
         time=infile.variables['time'][:]
         lats=infile.variables['lat'][:]
         lons=infile.variables['lon'][:]
         
         lats_orig = lats
         lons_orig = lons
         
         if args.roi:
            print("Extracting region of interest: ",min_lat, max_lat, min_lon, max_lon) 
            ext_lats = np.logical_and(lats >= float(min_lat), lats <= float(max_lat))
            lat_indices = np.where(ext_lats)
            lats = lats[lat_indices]
            ext_lons = np.logical_and(lons >= float(min_lon), lons <= float(max_lon))
            lon_indices = np.where(ext_lons)
            lons = lons[lon_indices]
            
            lat_len = lat_indices[0].shape
            lat_range = lat_indices[0][0:lat_len[0]-1]
            lon_len = lon_indices[0].shape
            lon_range = lon_indices[0][0:lon_len[0]-1]
            
            kd_blue_ext = kd_blue[0][lat_range[0]:lat_range[len(lat_range)-1],lon_range[0]:lon_range[len(lon_range)-1]]
            kd_blue_ext = kd_blue_ext[newaxis,...]
            kd_blue_shape = kd_blue_ext.shape
            print("Shape of Kd blue")
            print(kd_blue_shape)

            kd_green_ext = kd_green[0][lat_range[0]:lat_range[len(lat_range)-1],lon_range[0]:lon_range[len(lon_range)-1]]
            kd_green_ext = kd_green_ext[newaxis,...]
            kd_green_shape = kd_green_ext.shape
            print("Shape of Kd green")
            print(kd_green_shape)

            kd_red_ext = kd_red[0][lat_range[0]:lat_range[len(lat_range)-1],lon_range[0]:lon_range[len(lon_range)-1]]
            kd_red_ext = kd_red_ext[newaxis,...]
            kd_red_shape = kd_red_ext.shape
            print("Shape of Kd red")
            print(kd_red_shape)
            x_kd = kd_red_shape[1]
            y_kd = kd_red_shape[2]
         
            if args.seabed:
               print(z_bathy.shape)
               ext_bathy_lats = np.logical_and(bathy_lats >= float(min_lat), bathy_lats <= float(max_lat))
               bathy_lat_indices = np.where(ext_bathy_lats)
               ext_bathy_lons = np.logical_and(bathy_lons >= float(min_lon), bathy_lons <= float(max_lon))
               bathy_lon_indices = np.where(ext_bathy_lons)
               lat_len = bathy_lat_indices[0].shape
               print(lat_len)
               lat_range = bathy_lat_indices[0][0:lat_len[0]-1]
               lon_len = bathy_lon_indices[0].shape
               print(lon_len)
               lon_range = bathy_lon_indices[0][0:lon_len[0]-1]
               z_bathy_ext = z_bathy[0][lat_range[0]:lat_range[len(lat_range)-1],lon_range[0]:lon_range[len(lon_range)-1]]
               z_bathy_shape = z_bathy_ext.shape
               print("Resized bathymetry file shape: ", z_bathy_shape)
            
            ext_falchi_lats = np.logical_and(falchi_lats_orig >= float(min_lat), falchi_lats_orig <= float(max_lat))
            lat_indices = np.where(ext_falchi_lats)
            ext_falchi_lons = np.logical_and(falchi_lons_orig >= float(min_lon), falchi_lons_orig <= float(max_lon))
            lon_indices = np.where(ext_falchi_lons)
            
            lat_len = lat_indices[0].shape
            lat_range = lat_indices[0][0:lat_len[0]-1]
            lon_len = lon_indices[0].shape
            lon_range = lon_indices[0][0:lon_len[0]-1]
            
            falchi_ext = falchi_map[lat_range[0]:lat_range[len(lat_range)-1],lon_range[0]:lon_range[len(lon_range)-1]]
            falchi_lats = falchi_lats_orig[lat_range[0]:lat_range[len(lat_range)-1]]
            falchi_lons = falchi_lons_orig[lon_range[0]:lon_range[len(lon_range)-1]]
            print("Shape of Falchi map")
            print(falchi_ext.shape)
            xy_falchi = falchi_ext.shape
            x_falchi = xy_falchi[0]
            y_falchi = xy_falchi[1]
                        
         else:
            # extract data to match that of Falchi map - only need to account for latitude
            ext_lats = np.logical_and(lats >= miny, lats <= maxy)
            indices = np.where(ext_lats)
            lats = lats[indices]
         
            kd_blue_ext = kd_blue[0][indices][:]
            kd_blue_ext = kd_blue_ext[newaxis,...]
            kd_blue_shape = kd_blue_ext.shape

            kd_green_ext = kd_green[0][indices][:]
            kd_green_ext = kd_green_ext[newaxis,...]
            kd_green_shape = kd_green_ext.shape

            kd_red_ext = kd_red[0][indices][:]
            kd_red_ext = kd_red_ext[newaxis,...]
            kd_red_shape = kd_red_ext.shape
         
            if args.seabed:
               z_bathy_ext = z_bathy[0][indices][:]
               z_bathy_ext = z_bathy_ext[newaxis,...]
               z_bathy_shape = z_bathy_ext.shape
               print("Resized bathymetry file shape: ", z_bathy_shape)

         start = timer()
         if args.roi:
            print("Resizing IOP maps to ~1km - Falchi")
            #zoom_factor_lat = float(height)/float(kd_blue_shape[1])
            #zoom_factor_lon = float(width)/float(kd_blue_shape[2])
            zoom_factor_lat = float(y_falchi)/float(y_kd)
            zoom_factor_lon = float(x_falchi)/float(x_kd)
            print("Resizing Kd blue and changing fill value to -0.01")
            kd_blue_ext = np.ma.filled(kd_blue_ext.astype(float), -0.01)
            kd_blue_resize = zoom(kd_blue_ext, (1.0, zoom_factor_lat, zoom_factor_lon))
            kd_blue_ext = kd_blue_resize
            kd_blue_shape = kd_blue_resize.shape # This variable is needed with writing out the netCDF
            print("Resizing Kd green and changing fill value to -0.01")
            kd_green_ext = np.ma.filled(kd_green_ext.astype(float), -0.01)
            kd_green_resize = zoom(kd_green_ext, (1.0, zoom_factor_lat, zoom_factor_lon))
            kd_green_ext = kd_green_resize
            print("Resizing Kd red and changing fill value to -0.01")
            kd_red_ext = np.ma.filled(kd_red_ext.astype(float), -0.01)
            kd_red_resize = zoom(kd_red_ext, (1.0, zoom_factor_lat, zoom_factor_lon))
            kd_red_ext = kd_red_resize
            
            # no need to resize bathymetry - as this is done at the beginning by reading in a higher
            # resolution map
            #print("Resizing bathymetry matrix")
            
            # no need to resize falchi - but here for completeness
            falchi_resize = falchi_ext
            # resize of the latitude and longitude matrices
            lats = falchi_lats
            lons = falchi_lons
         else:
            # resize the Falchi map to 4km
            print("Resizing Falchi map to 4 km")
            zoom_factor_lat = float(kd_blue_shape[1])/float(height)
            zoom_factor_lon = float(kd_blue_shape[2])/float(width)
            falchi_resize = zoom(falchi_map, (zoom_factor_lat, zoom_factor_lon))

         end = timer()
         print("Time taken: (minutes)")
         print((end - start)/60.)
         
         falchi_resize = falchi_resize[newaxis,...]
         print("Shape of Falchi map")
         print(falchi_resize.shape)
         
         # calculate the above water values of irradiance in the blue, green, red
         # Falchi units are in mCd/m2
         # irradiance is in uW/m2
         print("Calculating surface irradiances")
         ##log_sfce_irr_blue_uW_m2 = m_blue*falchi_resize + c_blue
         ##log_sfce_irr_green_uW_m2 = m_green*falchi_resize + c_green
         ##log_sfce_irr_red_uW_m2 = m_red*falchi_resize + c_red
         
         ##sfce_irr_blue_uW_m2 = np.zeros(log_sfce_irr_blue_uW_m2.shape)
         ##sfce_irr_green_uW_m2 = np.zeros(log_sfce_irr_green_uW_m2.shape)
         ##sfce_irr_red_uW_m2 = np.zeros(log_sfce_irr_red_uW_m2.shape)
         ##realistic_blue = log_sfce_irr_blue_uW_m2 < 7.0
         ##realistic_green = log_sfce_irr_green_uW_m2 < 7.0
         ##realistic_red = log_sfce_irr_red_uW_m2 < 7.0

         ##sfce_irr_blue_uW_m2[realistic_blue] = 10**(log_sfce_irr_blue_uW_m2[realistic_blue])
         ##sfce_irr_green_uW_m2[realistic_green] = 10**(log_sfce_irr_green_uW_m2[realistic_green])
         ##sfce_irr_red_uW_m2[realistic_red] = 10**(log_sfce_irr_red_uW_m2[realistic_red])

         # Calculation with offset added (7/3/21).  This corroborated by 
         # working with the Tamir data in Eilat
         sfce_irr_blue_uW_m2 = m_blue*falchi_resize + c_blue
         sfce_irr_green_uW_m2 = m_green*falchi_resize + c_green
         sfce_irr_red_uW_m2 = m_red*falchi_resize + c_red
         
         # set surface values to zero where there is negative or zero signal from Falchi maps
         no_ALAN_blue = np.where(sfce_irr_blue_uW_m2 <= (c_blue+c_blue*0.01))
         no_ALAN_green = np.where(sfce_irr_green_uW_m2 <= (c_green+c_green*0.01))
         no_ALAN_red = np.where(sfce_irr_red_uW_m2 <= (c_red+c_red*0.01))

         sfce_irr_blue_uW_m2[no_ALAN_blue] = 0.
         sfce_irr_green_uW_m2[no_ALAN_green] = 0.
         sfce_irr_red_uW_m2[no_ALAN_red] = 0.
         
         print("Calculating irradiance ratios")
         sfce_irr_total_uW_m2 = sfce_irr_blue_uW_m2 + sfce_irr_green_uW_m2 + sfce_irr_red_uW_m2
         R_blue = np.zeros(sfce_irr_total_uW_m2.shape)
         R_green = np.zeros(sfce_irr_total_uW_m2.shape)
         R_red = np.zeros(sfce_irr_total_uW_m2.shape)
         
         R_non_zero = sfce_irr_total_uW_m2 > 0.
         
         R_blue[R_non_zero] = sfce_irr_blue_uW_m2[R_non_zero]/sfce_irr_total_uW_m2[R_non_zero]
         R_green[R_non_zero] = sfce_irr_green_uW_m2[R_non_zero]/sfce_irr_total_uW_m2[R_non_zero]
         R_red[R_non_zero] = sfce_irr_red_uW_m2[R_non_zero]/sfce_irr_total_uW_m2[R_non_zero]
         
         # calculate the total Kd (broadband) at 1 m depth ...
         print("Calculating total Kd")
         Kd_log_expression = np.ones(kd_blue_ext.shape)*np.nan
         Kd_total = np.ones(kd_blue_ext.shape)*np.nan
         Kd_log_expression = R_blue*np.exp(-1.0*kd_blue_ext) + R_green*np.exp(-1.0*kd_green_ext) + R_red*np.exp(-1.0*kd_red_ext)

         Kd_pos = np.logical_and(Kd_log_expression > 0., ~np.isnan(Kd_log_expression))
         Kd_total[Kd_pos] = -1.0*np.log(Kd_log_expression[Kd_pos])
         
         # mask out values which are at or below intercept
         #blue_mask = sfce_irr_blue_uW_m2 <= c_blue
         #green_mask = sfce_irr_green_uW_m2 <= c_green
         #red_mask = sfce_irr_red_uW_m2 <= c_red
         blue_mask = sfce_irr_blue_uW_m2 <= 0.
         green_mask = sfce_irr_green_uW_m2 <= 0.
         red_mask = sfce_irr_red_uW_m2 <= 0.
         
         sfce_irr_blue_uW_m2[blue_mask] = np.nan
         sfce_irr_green_uW_m2[green_mask] = np.nan
         sfce_irr_red_uW_m2[red_mask] = np.nan

         ind_mask = np.logical_and(blue_mask, green_mask, red_mask)

         # calculate the light level at the seafloor
         if args.seabed:
            print("Calculating irradiance at seafloor")
            irr_total_uW_m2 = np.ones(sfce_irr_total_uW_m2.shape)*np.nan
            if args.roi:
               z_bathy_ext = z_bathy_ext[newaxis,...]
            Kd_non_zero = Kd_total > 0.
            z_bathy_pos = z_bathy_ext > 0.
            sfce_irr_total_pos = sfce_irr_total_uW_m2 > 0.
            calc_mask = np.logical_and(Kd_non_zero, z_bathy_pos, sfce_irr_total_pos)
            #irr_total_uW_m2[Kd_non_zero] = sfce_irr_total_uW_m2[Kd_non_zero]*np.exp(-1.0*Kd_total[Kd_non_zero]*z_bathy_ext[Kd_non_zero])
            irr_total_uW_m2[calc_mask] = sfce_irr_total_uW_m2[calc_mask]*np.exp(-1.0*Kd_total[calc_mask]*z_bathy_ext[calc_mask])
            irr_total_uW_m2[ind_mask] = np.nan
         
         # calculate depth below which light levels drop below threshold
         if args.Zc:
            print("Calculating threshold depth")
            z_thresh = np.ones(sfce_irr_total_uW_m2.shape)*np.nan
            Kd_non_zero = np.logical_and(Kd_total > 0., ~np.isnan(Kd_log_expression))
            z_thresh[Kd_non_zero] = (-1.0/Kd_total[Kd_non_zero])*np.log(thresh_irr_total_uW_m2/sfce_irr_total_uW_m2[Kd_non_zero])
            z_thresh[ind_mask] = np.nan
         
         # output the resultant netCDF file containing Kd(RGB) and lat, lon
         w_nc_fid = nc.Dataset(outdir+ofname, 'w', format='NETCDF4')
         w_nc_fid.description = "Resized onto Falchi map dimensions spectral Kd values determined using IOP and radiative transfer"

         # Create the new dimensions from the resized variables
         data = {}
         new_dims = [kd_blue_shape[0], kd_blue_shape[1], kd_blue_shape[2]]
         counter = 0
         for dim in nc_dims:
            #w_nc_fid.createDimension(dim, infile.variables[dim].size)
            w_nc_fid.createDimension(dim, new_dims[counter])
            data[dim] = w_nc_fid.createVariable(dim, infile.variables[dim].dtype,(dim,))
            # You can do this step yourself but someone else did the work for us.
            for ncattr in infile.variables[dim].ncattrs():
               data[dim].setncattr(ncattr, infile.variables[dim].getncattr(ncattr))

            counter = counter+1

         # Assign the dimension data to the new NetCDF file.
         w_nc_fid.variables['time'][:] = time
         w_nc_fid.variables['lat'][:] = lats
         w_nc_fid.variables['lon'][:] = lons

         # Create output variables
         w_nc_var = w_nc_fid.createVariable('falchi_resize', 'f8', ('time', 'lat', 'lon'))
         w_nc_var.setncatts({'long_name': u"Falchi ALAN map",\
                       'units': u"mCd/m2", 'level_desc': u'Surface',\
                       'var_desc': u"Surface ALAN: mCd/m2"})
         w_nc_fid.variables['falchi_resize'][:] = falchi_resize
         
         # Create output variables
         #w_nc_var = w_nc_fid.createVariable('z_bathy_ext', 'f8', ('time', 'lat', 'lon'))
         #w_nc_var.setncatts({'long_name': u"Bathymetry",\
         #              'units': u"m", 'level_desc': u'Depth',\
         #              'var_desc': u"Bathymetry: m"})
         #w_nc_fid.variables['z_bathy_ext'][:] = z_bathy_ext 
         
         if args.Zc:
            w_nc_var = w_nc_fid.createVariable('z_thresh', 'f8', ('time', 'lat', 'lon'))
            w_nc_var.setncatts({'long_name': u"Depth below which light threshold reached",\
                          'units': u"m", 'level_desc': u'Surface',\
                          'var_desc': u"Threshold depth: m"})
            w_nc_fid.variables['z_thresh'][:] = z_thresh

            # Quick check to see if the Kd total values are being correctly calculated
            w_nc_var = w_nc_fid.createVariable('Kd_total', 'f8', ('time', 'lat', 'lon'))
            w_nc_var.setncatts({'long_name': u"Diffuse attenuation coefficient (Total)",\
                          'units': u"1/m", 'level_desc': u'Surface',\
                          'var_desc': u"Light Decay rate: 1/m"})
            w_nc_fid.variables['Kd_total'][:] = Kd_total

         if args.seabed:
            w_nc_var = w_nc_fid.createVariable('irr_total_uW_m2', 'f8', ('time', 'lat', 'lon'))
            w_nc_var.setncatts({'long_name': u"Total broadband irradiance at seafloor",\
                          'units': u"uW/m2", 'level_desc': u'Sea Floor',\
                          'var_desc': u"Irradiance: uW/m2"})
            w_nc_fid.variables['irr_total_uW_m2'][:] = irr_total_uW_m2
         

         w_nc_fid.close()  # close the new file
         infile.close()
   sys.exit()
      
if __name__=='__main__':
   main()
