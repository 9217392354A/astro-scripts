import shutil
import os
from os.path import join as pj

folder = "/home/fossil/spx6cff/coma/"
mvFile = "Hipe_rename_files.txt"
cpFile = "Hipe_move_files.txt"
zpFile = "Hipe_unzip_files.txt"
allFile = "Hipe_rename_cp_unzip.txt"

readDir = "/home/fossil/spx6cff/coma/hipedata/cfuller132117751/"
writeDirBlue = "/home/fossil/spx6cff/coma/hipedata/map/PACSBLUE/"
writeDirRed = "/home/fossil/spx6cff/coma/hipedata/map/PACSRED/"

topDirs = os.listdir(readDir)
topDirs = [x for x in topDirs if x[-4:] != ".xml"]

outFile = open(pj(folder,mvFile), 'w')
outFile2 = open(pj(folder,cpFile), 'w')
outFile3 = open(pj(folder,zpFile), 'w')
outFile4 = open(pj(folder,allFile), 'w')

count  = 0
#blue
for currentDir in topDirs:
    count += 1
    
    path = str(readDir) + str(currentDir) + "/level2/HPPPMAPB/"
    fileName = str(os.listdir(path))
    fileName = fileName.split("'")
    fileName = str(fileName[1])
    destPath = str(writeDirBlue)
    fileNameNew = "blue" + str(count) + ".fits.gz"
    line = "mv -v " + path + str(fileName) +  " " + path + str(fileNameNew)+ "\n"
    line2 = "cp -v " + path + str(fileName) +  " " + destPath + str(fileNameNew) + "\n"
    line3 = "gunzip -v " + destPath +  str(fileNameNew) + "\n"
    outFile.write(line)
    outFile2.write(line2)
    outFile3.write(line3)
    outFile4.write(line + line2 +line3)
    
    
count2 = 0
#red
for currentDir in topDirs:
    count2 += 1
    
    path = str(readDir) + str(currentDir) + "/level2/HPPPMAPR/"
    fileName = str(os.listdir(path))
    fileName = fileName.split("'")
    fileName = str(fileName[1])
    destPath = str(writeDirRed)
    fileNameNew = "red" + str(count2) + ".fits.gz"
    line = "mv -v " + path + str(fileName) +  " " + path + str(fileNameNew)+ "\n"
    line2 = "cp -v " + path + str(fileName) +  " " + destPath + str(fileNameNew) + "\n"
    line3 = "gunzip -v " + destPath +  str(fileNameNew)+ "\n"
    outFile.write(line)
    outFile2.write(line2)
    outFile3.write(line3)
    outFile4.write(line + line2 +line3)

print "examples"
print line
print line2
print line3
print "Number of files...Blue:", count, " Red:", count2 

outFile.close()
outFile2.close()
outFile3.close()
outFile4.close()