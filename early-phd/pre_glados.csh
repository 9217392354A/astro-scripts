#! /bin/csh -f
#Chris Fuller Jan 2011
#Program that can be used to prep cubes for use with source finder GLADOS(Rhys Talyor)
#
#instructions
#
#edit imputs to your cube names and then run
#
#make sure that you are running this in the same directory as your files are located and that you first type "miriadstart"
#

#inputs ...

set pola = northpola_withmasking_poly2.fits
set polb = northpolb_withmasking_poly2.fits
set main = north_withmasking_poly2.fits

##################################################################################

echo Stage 1 converting into miriad format 
#Converting to .mir format to keep miriad sweet

fits in=$pola op=xyin out=pola.mir
fits in=$polb op=xyin out=polb.mir
fits in=$main op=xyin out=main.mir

echo Stage 2 performing reorder
#Reorder axis to make hanning faster

reorder in=pola.mir mode=321 out=polare.mir
reorder in=polb.mir mode=321 out=polbre.mir
reorder in=main.mir mode=321 out=mainre.mir

echo Stage 3 perfoming hanning smoothing to all cubes with a width of 3
#Smoothing all cubes by 3

hanning in=polare.mir width=3 out=polareh3.mir
hanning in=polbre.mir width=3 out=polbreh3.mir
hanning in=mainre.mir width=3 out=mainreh3.mir

echo Stage 4 further hanning smoothing to width 5 to both porlarisations
#smoothing pol cubes by further 5

hanning in=polareh3.mir width=5 out=polareh35.mir
hanning in=polbreh3.mir width=5 out=polbreh35.mir

echo Stage 5 reording back
#reordering back to original axis

reorder in=polareh35.mir mode=321 out=polah35.mir
reorder in=polbreh35.mir mode=321 out=polbh35.mir
reorder in=mainreh3.mir mode=321 out=mainh3.mir

echo converting back to fits format
#Converting back to .fits format to keep GLADOS happy

fits in=polah35.mir op=xyout out=GLADOS_$pola
fits in=polbh35.mir op=xyout out=GLADOS_$polb
fits in=mainh3.mir op=xyout out=GLADOS_$main

echo Program finished
echo outputs prefixed with GLADOS all .mir files can be deleted....
