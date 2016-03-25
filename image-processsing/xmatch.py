#Program created by Chris Fuller May 2012
#Program to cross match two catalogues

# import stuff
from numpy import *
import numpy
import scipy
import math
import sys
import os
from os.path import join as pj

#User varibles
threshold = 10.0 / 3600.0
#Col where ra appears in file for each respective table
Acol = 1
Bcol = 2

#File stuff
folder = "/Users/chrisfuller/Dropbox/coma/Catalogues/"
fileA = "comasdssroughSLIM.csv"
fileB = "goldmineSLIM.csv"
outfile = "matched"

#inport files into arrays with headers
firstFile = loadtxt(pj(folder,fileA), dtype=str, skiprows=0, delimiter=",", unpack=False)
secondFile = loadtxt(pj(folder,fileB), dtype=str, skiprows=0, delimiter=",", unpack=False)
headA = firstFile[0,:]
headB = secondFile[0,:]
catA = firstFile[1:len(firstFile),:]
catB = secondFile[1:len(secondFile),:]
#Shape array so that it has a single row that can be added too
head = append(headA,headB)
head = head.reshape(1,len(head))
##################################################################################
#functions 

#xmatch is a function to look through arrayA for matchs and then adding thsese to the new table
#however where there are no matchs adding zero's 
def xmatch(arrayA,arrayB,Acol,Bcol,head):  
    #first looping through source cat and then looping through for each source to test against 
    for rowA in arrayA:
        count = 0
        selection = []
        for rowB in arrayB:
            raA = float(rowA[Acol])
            decA = float(rowA[(Acol + 1)])
            raB = float(rowB[Bcol])
            decB = float(rowB[(Bcol+ 1)])
            #see if source is within threshold included ing cos dec projection term            
            d = sqrt((abs(raA - raB)*math.cos(decA))**2.0 + abs((abs(decA - decB)*math.cos(decB))**2.0 ))
            if d < threshold:
                count += 1
                selection = rowB
                print selection   
        if count == 0:
            raise ' need to add zeros here !'
            #add 0 to output array that is the same len as rowB
        elif count == 1:
            #wrting selection to array and then making then adding row's together then adding a new row
            #to the main table
            selection = append(rowA,selection)
            selection = selection.reshape(1,len(selection)) 
            head = append(head,selection, axis=0)
        else:
            raise 'two sources found within threshold'       
    return head
            
##################################################################################

matches = xmatch(catA,catB,Acol,Bcol,head)
   
    