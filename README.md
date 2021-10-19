# ALAN
Artificial Light at Night

Code to calculate the critical depth (based on the light sensitivity of Calanus) of artificial light penetration in the global ocean.

Methodology outlined in Smyth et al., (2021) "A global atlas of artificial light at night under the sea" 

Python code: calculate_ALAN.py 
------------------------------
Calculates the Zc ALAN for a given region of interest (ROI) defined as min lat, max lat, min lon, max lon.  Outputs as netCDF file


Shell script: batch_regional_tiles.csh
--------------------------------------
Creates a 10 x 10 degree tile if there is coastline within the tile (i.e. doesn't create anything for the middle of the Pacific Ocean!)

Separate routine required if need to create a composite of several tiles (e.g. an entire continental margin) 


