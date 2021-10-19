#!/bin/csh 
#ftp://wcobuoy:Toogh0shaephai@ftp.rsg.pml.ac.uk/
setenv DISPLAY :3
set exdir = /users/rsg/tjsm/ALAN

#### UK ####
set maxlimit=30      #| in plot_ncfiles.py
set irrmaxlimit=2.0  
#nice +19 $exdir/calculate_ALAN.py --Zc --roi 45 65 -10 10
nice +19 $exdir/plot_ncfiles.py --Zc --roi 45 65 -10 10 --idir /users/rsg/tjsm/ALAN/Global/SciAdv --odir /users/rsg/tjsm/ALAN/Global/SciAdv/ 

#### Persian Gulf ###
#set maxlimit=50      # in plot_ncfiles.py
#set irrmaxlimit=5.0  #
#nice +19 $exdir/calculate_ALAN.py --Zc --roi 22 32 48 58
#nice +19 $exdir/plot_ncfiles.py --Zc --roi 22 32 48 58 --idir /users/rsg/tjsm/tmp/Zc --odir /users/rsg/tjsm/tmp/Zc/ 

#montage -geometry 400x800 ESACCI-OC-MAPPED-CLIMATOLOGY-1M_MONTHLY_4km_GEO_PML_OCx_QAA-Zcritical-clear-05-fv4.0_22S_32N_48W_58E.png ESACCI-OC-MAPPED-CLIMATOLOGY-1M_MONTHLY_4km_GEO_PML_OCx_QAA-Zcritical-clear-05-fv4.0_45S_65N_-10W_10E.png Smyth_figure2.png
#montage -geometry 400x800 ESACCI-OC-MAPPED-CLIMATOLOGY-1M_MONTHLY_4km_GEO_PML_OCx_QAA-Zcritical-clear-05-fv4.0_22S_32N_48W_58E.png ESACCI-OC-MAPPED-CLIMATOLOGY-1M_MONTHLY_4km_GEO_PML_OCx_QAA-Zcritical-clear-05-fv4.0_45S_65N_-10W_10E.png Smyth_figure2.png
exit(0)
