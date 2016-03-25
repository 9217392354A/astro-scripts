#This program read in the positions of every H-ATLAS source, which has
#already been matched to a VLA source, an also reads in the positions
#of every SDF source. It then performs a simple enarest neighbour routine
#to find the best matches for every H-ATLAS source. This replaces the
#likelihood method for matching these catalogues which is not right to
#use.

#Matthew Allen - Cardiff University
#Sep 2013

####################### Import Stuff #######################
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import join as pj

####################### Parameters #######################

pi = 2.0*(np.arccos(0.0))

####################### Variables #######################
MaxSearchRadius = 1.0
RandomPositionsLoop = False

####################### Data #################################
folder_Primary = "/home/gould-belt/c0716754/Work/NGPGalaxyMSWork/EzzysRedshiftEstimator/" #H-ATLAS sources
folder_Secondary = "/home/gould-belt/c0716754/Work/Data/SDF/" #SDF sources

datafile_Primary = "NGP_LikRelZ.dat"
datafile_Secondary = "SDF_FullPhoto.dat"

SourceDataPrimary = np.loadtxt( pj(folder_Primary,datafile_Primary), comments='#', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19))
SourceDataPrimaryIDs = np.genfromtxt( pj(folder_Primary,datafile_Primary), dtype='str', comments='#', usecols=(20))
SourceDataSecondary = np.loadtxt( pj(folder_Secondary,datafile_Secondary), comments='#')

folder_Plots = "/home/gould-belt/c0716754/Work/NGPGalaxyMSWork/RedshiftEstimation/Plots"


FileOut = open('ReliableSources.dat','w')
FileOutForHyperZ = open('UBVRI_SDF.cat','w')
FileHeader = "# HerschelRA HerschelDec RadioRA RadioDec SDFRA SDFDec CounterpartRadius HerschelID PrimaryLoop f250/e250 f250 f350 f500 e250 e350 e500 EzzyZ BMag Berr VMag VErr RMag RErr IMag Ierr ZMag ZErr JMag Jerr KMag Kerr photoZaper photoZauto RadioFlux RadioFluxError Reliability SDFOrNot SDFHerschelCounterpartRadius H-ATLASID" + "\n"
			
FileOut.write(FileHeader)


#---- Now create a new array, but only with the primary sources that are covered by the SDF map.

SDFSelection = np.where((SourceDataPrimary[:,PrimaryRA] > 200.8816) & (SourceDataPrimary[:,PrimaryRA] < 201.444) & (SourceDataPrimary[:,PrimaryDec] > 27.182) & (SourceDataPrimary[:,PrimaryDec] < 27.799))

#SDFSelection = np.where((SourceDataPrimary[:,PrimaryRA] > 200.8816))

CutSourceDataPrimary = SourceDataPrimary[SDFSelection[0],:]

CutSourceDataPrimaryIDs = SourceDataPrimaryIDs[SDFSelection[0],:]

print CutSourceDataPrimaryIDs.shape[0]


#---- Find the lengths of the Primary and Secondary arrays
NumberOfPrimarySources = SourceDataPrimary.shape[0]
NumberOfCutPrimarySources = CutSourceDataPrimary.shape[0]
print 'No. Prim sources: ' + str(NumberOfPrimarySources)
print 'No. cut Prim sources: ' + str(NumberOfCutPrimarySources)

NumberOfSecondarySources = SourceDataSecondary.shape[0]
print 'No. Sec sources: ' + str(NumberOfSecondarySources)

NumberOfRandomPositions = 1000

#---- Generate an array of random positions
RandomPositionsArray = np.random.rand(NumberOfRandomPositions,2)

for i in range (0,NumberOfRandomPositions):
	RandomPositionsArray[i,0] = (RandomPositionsArray[i,0]*0.5624) + 200.8816
	
for i in range (0,NumberOfRandomPositions):
	RandomPositionsArray[i,1] = (RandomPositionsArray[i,1]*0.617) + 27.1820
	

#---- Loop over Primary and Secondary sources
#NumberOfPrimarySources = 10

for PrimaryLoop in range(0,NumberOfPrimarySources):
	PrimaryNearbysourceCounter = 0
	
	for SecondaryLoop in range(0,NumberOfSecondarySources):
	
		DeltaRA = CutSourceDataPrimary[PrimaryLoop,PrimaryRA] - SourceDataSecondary[SecondaryLoop,SecondaryRA]
		DeltaDec = CutSourceDataPrimary[PrimaryLoop,PrimaryDec] - SourceDataSecondary[SecondaryLoop,SecondaryDec]
		
		CounterpartRadius = np.sqrt( ((DeltaRA*np.cos(SourceDataSecondary[SecondaryLoop,SecondaryDec]* (pi/180.0)))**2.0 ) + (DeltaDec**2.0) )		
		CounterpartRadius = CounterpartRadius * 3600.0
		
#		print CounterpartRadius

		if (CounterpartRadius < MaxSearchRadius):
			PrimaryNearbysourceCounter += 1

#---- General output file
#---- Output: HerschelRA, HerschelDec, RadioRA, RadioDec, SDFRA, SDFDec, CounterpartRadius, HerschelID, PrimaryLoop, f250/e250, f250, f350, f500, e250, e350, e500, EzzyZ, BMag, Berr, VMag, VErr, RMag, RErr, IMag, Ierr, ZMag, ZErr, JMag, Jerr, KMag, Kerr, photoZ aper, photoZ auto, RadioFlux,

			lineReliableSources = '{0:7.3f} {1:7.3f} {2:7.3f} {3:7.3f} {4:7.3f} {5:7.3f} {6:7.3f} {7:4d} {8:4d} {9:7.3f} {10:7.4f} {11:7.4f} {12:7.4f} {13:7.4f} {14:7.4f} {15:7.4f} {16:7.3f} {17:7.3f} {18:7.3f} {19:7.3f} {20:7.3f} {21:7.3f} {22:7.3f} {23:7.3f} {24:7.3f} {25:7.3f} {26:7.3f} {27:7.3f} {28:7.3f} {29:7.3f} {30:7.3f} {31:7.3f} {32:7.3f} {33:7.3f} {34:7.3f} {35:8.8s} {36:2d} {37:7.3f} {38:7.24s} {39}'.format((CutSourceDataPrimary[PrimaryLoop,1]), (CutSourceDataPrimary[PrimaryLoop,2]), (CutSourceDataPrimary[PrimaryLoop,PrimaryRA]), (CutSourceDataPrimary[PrimaryLoop,PrimaryDec]), (SourceDataSecondary[SecondaryLoop,SecondaryRA]), (SourceDataSecondary[SecondaryLoop,SecondaryDec]), (CutSourceDataPrimary[PrimaryLoop,12]), (int(CutSourceDataPrimary[PrimaryLoop,0])), int(PrimaryLoop), (CutSourceDataPrimary[PrimaryLoop,5]/CutSourceDataPrimary[PrimaryLoop,9]), (CutSourceDataPrimary[PrimaryLoop,5]), (CutSourceDataPrimary[PrimaryLoop,6]), (CutSourceDataPrimary[PrimaryLoop,7]), (CutSourceDataPrimary[PrimaryLoop,8]), (CutSourceDataPrimary[PrimaryLoop,9]), (CutSourceDataPrimary[PrimaryLoop,10]), (CutSourceDataPrimary[PrimaryLoop,18]), (SourceDataSecondary[SecondaryLoop,2]), (SourceDataSecondary[SecondaryLoop,3]), (SourceDataSecondary[SecondaryLoop,4]), (SourceDataSecondary[SecondaryLoop,5]), (SourceDataSecondary[SecondaryLoop,6]), (SourceDataSecondary[SecondaryLoop,7]), (SourceDataSecondary[SecondaryLoop,8]), (SourceDataSecondary[SecondaryLoop,9]), (SourceDataSecondary[SecondaryLoop,10]), (SourceDataSecondary[SecondaryLoop,11]), (SourceDataSecondary[SecondaryLoop,12]), (SourceDataSecondary[SecondaryLoop,13]), (SourceDataSecondary[SecondaryLoop,14]), (SourceDataSecondary[SecondaryLoop,15]), (SourceDataSecondary[SecondaryLoop,16]), (SourceDataSecondary[SecondaryLoop,17]), (CutSourceDataPrimary[PrimaryLoop,16]), (CutSourceDataPrimary[PrimaryLoop,17]), str(CutSourceDataPrimary[PrimaryLoop,15]), 1, CounterpartRadius, (CutSourceDataPrimaryIDs[PrimaryLoop]), "\n")		
	
			

			FileOut.write(lineReliableSources)
			
#---- Output for HyperZ			
#---- Output: HerschelIDBMag,VMag, RMag, IMag, ZMag, JMag, KMag, Berr, VErr, RErr, Ierr, ZErr, Jerr, Kerr
		
#			lineHyperZ = str(int(CutSourceDataPrimary[PrimaryLoop,0])) + "  " + str(SourceDataSecondary[SecondaryLoop,2]) + "  " + str(SourceDataSecondary[SecondaryLoop,4]) + "  " + str(SourceDataSecondary[SecondaryLoop,6]) + "  " + str(SourceDataSecondary[SecondaryLoop,8]) + "  " + str(SourceDataSecondary[SecondaryLoop,10]) + "  " + str(SourceDataSecondary[SecondaryLoop,12]) + "  " + str(SourceDataSecondary[SecondaryLoop,14]) + "  " + str(SourceDataSecondary[SecondaryLoop,3]) + "  " + str(SourceDataSecondary[SecondaryLoop,5]) + "  " + str(SourceDataSecondary[SecondaryLoop,7]) + "  " + str(SourceDataSecondary[SecondaryLoop,9]) + "  " + str(SourceDataSecondary[SecondaryLoop,11]) + "  " + str(SourceDataSecondary[SecondaryLoop,13]) + "  " + str(SourceDataSecondary[SecondaryLoop,15]) + "\n"
			
#			FileOutForHyperZ.write(lineHyperZ)
		
	if(PrimaryNearbysourceCounter > 1):
		print CutSourceDataPrimary[PrimaryLoop,0], PrimaryNearbysourceCounter
			
			
#---- Print out data for H-ATLAS sources without a SDF match			
	if(PrimaryNearbysourceCounter == 0):
	

		lineReliableSources = '{0:7.3f} {1:7.3f} {2:7.3f} {3:7.3f} {4:7.3f} {5:7.3f} {6:7.3f} {7:4d} {8:4d} {9:7.3f} {10:7.4f} {11:7.4f} {12:7.4f} {13:7.4f} {14:7.4f} {15:7.4f} {16:7.3f} {17:7.3f} {18:7.3f} {19:7.3f} {20:7.3f} {21:7.3f} {22:7.3f} {23:7.3f} {24:7.3f} {25:7.3f} {26:7.3f} {27:7.3f} {28:7.3f} {29:7.3f} {30:7.3f} {31:7.3f} {32:7.3f} {33:7.3f} {34:7.3f} {35:8.8s} {36:1d} {37:7.3f} {38:7.24s} {39}'.format((CutSourceDataPrimary[PrimaryLoop,1]), (CutSourceDataPrimary[PrimaryLoop,2]), (CutSourceDataPrimary[PrimaryLoop,PrimaryRA]), (CutSourceDataPrimary[PrimaryLoop,PrimaryDec]), (SourceDataSecondary[SecondaryLoop,SecondaryRA]), (SourceDataSecondary[SecondaryLoop,SecondaryDec]), (CutSourceDataPrimary[PrimaryLoop,12]), (int(CutSourceDataPrimary[PrimaryLoop,0])), int(PrimaryLoop), (CutSourceDataPrimary[PrimaryLoop,5]/CutSourceDataPrimary[PrimaryLoop,9]), (CutSourceDataPrimary[PrimaryLoop,5]), (CutSourceDataPrimary[PrimaryLoop,6]), (CutSourceDataPrimary[PrimaryLoop,7]), (CutSourceDataPrimary[PrimaryLoop,8]), (CutSourceDataPrimary[PrimaryLoop,9]), (CutSourceDataPrimary[PrimaryLoop,10]), (CutSourceDataPrimary[PrimaryLoop,18]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (CutSourceDataPrimary[PrimaryLoop,16]), (CutSourceDataPrimary[PrimaryLoop,17]), str(CutSourceDataPrimary[PrimaryLoop,15]), 0, -99.0, (CutSourceDataPrimaryIDs[PrimaryLoop]), "\n")					
						
		FileOut.write(lineReliableSources)				

#---- Loop over Random positions and Secondary Sources
#NumberOfPrimarySources = 10

if (RandomPositionsLoop == True):

	PrimaryNearbysourceCounter = 0

	for PrimaryLoop in range(0,NumberOfRandomPositions):


		for SecondaryLoop in range(0,NumberOfSecondarySources):
	
			DeltaRA = RandomPositionsArray[PrimaryLoop,0] - SourceDataSecondary[SecondaryLoop,SecondaryRA]
			DeltaDec = RandomPositionsArray[PrimaryLoop,1] - SourceDataSecondary[SecondaryLoop,SecondaryDec]
		
			CounterpartRadius = np.sqrt( ((DeltaRA*np.cos(SourceDataSecondary[SecondaryLoop,SecondaryDec]* (pi/180.0)))**2.0 ) + (DeltaDec**2.0) )		
			CounterpartRadius = CounterpartRadius * 3600.0

			if (CounterpartRadius < MaxSearchRadius):
				PrimaryNearbysourceCounter += 1
				print CounterpartRadius

	print "Number of nearby random sources: "		
	print PrimaryNearbysourceCounter
	print float(PrimaryNearbysourceCounter)/float(NumberOfRandomPositions)

		
print 'done'

FileOut.close()
		
