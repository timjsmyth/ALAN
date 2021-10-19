#!/bin/csh 
#ftp://wcobuoy:Toogh0shaephai@ftp.rsg.pml.ac.uk/
setenv DISPLAY :3
set exdir = /users/rsg/tjsm/ALAN

set start_lat = -50
set end_lat = 70
set start_lon = -180
set end_lon = 180

set lat = $start_lat
set ntiles = 1
set ntiles_to_process = 0

while ($lat < $end_lat)
   set slat = $lat 
   @ lat = $lat + 10
   set elat = $lat
   set lon = $start_lon
   while ($lon < $end_lon)
      set slon = $lon
      @ lon = $lon + 10
      set elon = $lon
      echo "==================="
      echo "Tile number: $ntiles (of 432)"
      echo $slat $elat $slon $elon
      set coast_flag = `$exdir/Global/mask/check_for_coast.py --roi $slat $elat $slon $elon`
      @ ntiles_to_process = $ntiles_to_process + $coast_flag
      if ($coast_flag) then
         nice +19 $exdir/calculate_ALAN.py --Zc --roi $slat $elat $slon $elon --odir /run/media/tjsm/Seagate\ Portable\ Drive/Seagate/data/ALAN/Zc/
         #nice +19 $exdir/plot_ncfiles.py --Zc --roi $slat $elat $slon $elon --idir /run/media/tjsm/Seagate\ Portable\ Drive/Seagate/data/ALAN/Zc/ --odir /run/media/tjsm/Seagate\ Portable\ Drive/Seagate/data/ALAN/Zc/png/
      endif
      echo "Total number processed: $ntiles_to_process"
      @ ntiles++
   end
end 
@ ntiles--

echo "Total number of tiles: $ntiles"

exit(0)
