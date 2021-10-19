#!/usr/bin/env python
#-*- coding:utf-8 -*-

import netCDF4
import numpy
import pylab
import os
import sys
from mpl_toolkits.basemap import Basemap
import argparse
import warnings
warnings.simplefilter('ignore', numpy.RankWarning)

def main():

   datadir='/home/scratch/data/workspace/ALAN/paper/Zc/'
   outdir='/home/scratch/data/workspace/ALAN/paper/Zc/png/' 

   parser = argparse.ArgumentParser(description='Produce ALAN pngs from netCDF')

   parser.add_argument("--Zc",  action='store_true',  default=False, help="Output critical depth pngs.")
   parser.add_argument("--seabed",  action='store_true',  default=False, help="Output irradiance at the seafloor pngs.")
   parser.add_argument("--roi", type=int, nargs=4, help="ROI - Region of interest")
   parser.add_argument("--idir", type=str, default='', help="Input directory for netcdf imagery")
   parser.add_argument("--odir", type=str, default='', help="Output directory for png imagery")
   parser.add_argument("--all", action='store_true', default=False, help="Process all files in directory.")
   parser.add_argument("--no_falchi", action='store_true', default=False,help="No Falchi map available.")
   args = parser.parse_args()

   pylab.rc('font', size=12)          # controls default text sizes
   #pylab.rc('axes', titlesize=12)     # fontsize of the axes title
   pylab.rc('axes', labelsize=18)    # fontsize of the x and y labels
   pylab.rc('xtick', labelsize=16)    # fontsize of the tick labels
   pylab.rc('ytick', labelsize=16)    # fontsize of the tick labels
   pylab.rc('legend', fontsize=16)    # legend fontsize
   pylab.rc('figure', titlesize=16)  # fontsize of the figure title

   if args.roi:
      print("Region of interest selected")
      roi_box = args.roi
      print("ROI: ", roi_box)
      min_lat = roi_box[0]
      max_lat = roi_box[1]
      min_lon = roi_box[2]
      max_lon = roi_box[3]
      filestr_box = '_'+str(min_lat)+'S_'+str(max_lat)+'N_'+str(min_lon)+'W_'+str(max_lon)+'E'
      print(filestr_box)

   # prevent both the seabed and critical depth options being selected
   if args.seabed and args.Zc:
      print("Select one of either --Zc or --seabed")
      sys.exit()

   # Critical depth
   if args.Zc:
      print("Creating critical depth pngs")
      if args.idir:
         datadir=args.idir
      else:
         datadir='/home/scratch/data/workspace/ALAN/paper/Zc/'
      if args.odir:
         outdir = args.odir
      else: 
         outdir=datadir+'png/' 
   # Irradiance at the seafloor
   if args.seabed:
      print("Creating seafloor irradiance pngs")
      if args.idir:
         datadir=args.idir
      else:
         datadir='/home/scratch/data/workspace/ALAN/seafloor/'
      if args.odir:
         outdir = args.odir
      else:
         outdir=datadir+'png/' 
   
   # if a roi is selected - want the coastline in high resolution
   # default is low
   resolution = 'l'
   if args.roi:
      print("High resolution coastline selected")
      resolution = 'f'
      
   # load in the data
   for filename in sorted(os.listdir(datadir)):
      print(filename)
      if args.roi:
         endfilestr = "fv4.0"+filestr_box+".nc"
      else:
         endfilestr = "fv4.0.nc"
      if args.all:
         endfilestr = ".nc"
      print(endfilestr)
      if filename.endswith(endfilestr): 
         fname = os.path.join(datadir, filename)
         print("Input filename: ", filename)
         nc=netCDF4.Dataset(fname)
         if args.Zc:
            zdata=nc.variables['z_thresh'][:][0]
         if args.seabed:
            zdata=nc.variables['irr_total_uW_m2'][:][0]
         if args.no_falchi:
            irrdata=zdata
         else:
            irrdata = nc.variables['falchi_resize'][:][0]
         lats=nc.variables['lat'][:]
         lons=nc.variables['lon'][:]
         nc.close()
         getext = filename.split('.')
         ofname = getext[0]+'.'+getext[1]+'.png'
         print(" Output filename: ", ofname)
         # get a stretch for colour scale - use 1% and 99% percentiles
         zflat=zdata.flatten() # flatten it to 1D
         irrflat = irrdata.flatten() # flatten it to 1D

         if args.seabed:
            zflat = numpy.log10(zflat)
            zdata = numpy.log10(zdata)
            
         I=numpy.where(~numpy.isnan(zflat)) # find which data are not NaN
         if len(I) > 100:
            minlimit=numpy.percentile(zflat[I],1) # get a lower scaling limit
            maxlimit=numpy.percentile(zflat[I],99) # upper scaling limit
         else:
            minlimit=0
            maxlimit=20
         zdata=numpy.flipud(zdata) # else it appears upside down in the basemap image

         I=numpy.where(irrdata > 0.) # find which data are not NaN

         if len(I) > 100:
            irrminlimit=numpy.percentile(irrdata[I],1) # get a lower scaling limit
            irrmaxlimit=numpy.percentile(irrdata[I],99) # upper scaling limit
         else:
            irrminlimit=0.
            irrmaxlimit=2.0
            
         irrdata=numpy.flipud(irrdata)
         
         if args.Zc:
            I = numpy.where(numpy.isnan(zdata))
            zdata[I] = 0.0
            minlimit=0.0
            maxlimit=30
            
         if args.seabed:
            minlimit = -2.0
            maxlimit = 1.0
            
         # create a basic image
         if args.roi:
            fig = pylab.figure(figsize=(10,20)) # change the size so that it is zoomable for small details
         else:
            fig = pylab.figure(figsize=(10,10)) # change the size so that it is zoomable for small details

         if args.no_falchi:
            fig1 = fig.add_subplot(1,1,1)
            a_map=Basemap(llcrnrlon=lons.min(),llcrnrlat=lats.min(),urcrnrlon=lons.max(),urcrnrlat=lats.max(),epsg=4326,resolution=resolution,suppress_ticks=False)
            a_map.drawcoastlines(color='1.0')
            a_map.fillcontinents(color='0.0')

            a_map.imshow(zdata,vmin=minlimit,vmax=maxlimit)
            pylab.colorbar(label='Critical depth [m]',orientation='horizontal') # add a colour bar
            
         else: 
            fig1 = fig.add_subplot(2,1,1)
            a_map=Basemap(llcrnrlon=lons.min(),llcrnrlat=lats.min(),urcrnrlon=lons.max(),urcrnrlat=lats.max(),epsg=4326,resolution=resolution,suppress_ticks=False)
            a_map.imshow(irrdata, vmin=irrminlimit, vmax=irrmaxlimit)
            a_map.drawcoastlines(color='1.0')
            pylab.colorbar(label='ALAN sky brightness (mcd m$^{-2}$)') # add a colour bar
            pylab.xlabel("Longitude ($^\circ$)")
            pylab.ylabel("Latitude ($^\circ$)")

            fig2 = fig.add_subplot(2,1,2)
            a_map=Basemap(llcrnrlon=lons.min(),llcrnrlat=lats.min(),urcrnrlon=lons.max(),urcrnrlat=lats.max(),epsg=4326,resolution=resolution,suppress_ticks=False)
            if args.roi: 
               a_map.drawcoastlines(color='1.0')
               a_map.fillcontinents(color='0.0')
            else: 
               a_map.drawcoastlines(color='0.0')

            a_map.imshow(zdata,vmin=minlimit,vmax=maxlimit)

            if args.Zc:
               pylab.colorbar(label='Critical depth (m)') # add a colour bar
            if args.seabed:
               pylab.colorbar(label='Log10(Irradiance) at seafloor [uW/m2]') # add a colour bar
            
         pylab.xlabel("Longitude ($^\circ$)")
         pylab.ylabel("Latitude ($^\circ$)")

         print("Saving figure")
         pylab.savefig(outdir+ofname,dpi=300) # save to 'publication' standard resolution
      #sys.exit()
   sys.exit()
      
if __name__=='__main__':
   main()


