#! /bin/csh -f
#
# script to concatenate fits cubes 
#
set i=1
set label = 00$i
while (-e $1_$label.fits)
  fits in=$1_$label.fits op=xyin out=$1_$label.mir
  @ i+=1
  if ( $i < 10 ) then
    set label = 00$i
  else if ( $i < 100 ) then
    set label = 0$i
  else
    set label = $i
  endif
end
if ( $i == 1 ) then
  echo $1_$label.fits not found
else
  imcat "in=$1_???.mir" out=$1.mir
  fits in=$1.mir op=xyout out=$1.fits
  echo Output miriad cube = $1.mir
  echo Output fits cube = $1.fits
endif
