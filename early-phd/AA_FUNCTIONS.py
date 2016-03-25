# This program is solely for storing functions and thier decriptions

##############################    Catalogue work ##################################
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

 ########################### input/output ##########################################
 
 # for saving an array to a text file .csv
def savetext(folder,filename,array):
    outfile = open(pj(folder, filename), 'w')
    for row in array:
        count = 0
        for element in row:
            count += 1
            if count < len(row):
                x = element + ","
            else:
                x = element
            outfile.write(x)
        outfile.write("\n")
    outfile.close()
    

############################     Fitting        #####################################   
# Function for genorating a exponential curve    
def genCurve(x,xDat):
    y =[]
    A = x[0]
    B = x[1]
    C = x[2]
    
    for i in range(0,len(xDat)):
        y.append(A*math.exp(-xDat[i]/B)+ C)
    return y
    
#function for fitting an exponential curve  
def curve(x,*args):    
    xDat = args[0]
    yDat = args[1]
    yError=args[2]
    A = x[0]
    B = x[1]
    C = x[2]
    ci = 0.0
    for i in range(0,len(xDat)):
        fx = A*math.exp(-xDat[i]/B)+ C
        ci += ((yDat[i] - fx)**2)/ (yError[i]**2)
    return ci
    
#Program for fitting a gaussian distribution 
def gauss(x,*args):    
    xDat = args[0]
    yDat = args[1]
    yError=args[2]
       
    hmax = x[0]
    hmin = x[1]
    mean =x[2]
    sigma = x[3]

    ci = 0.0
    for i in range(0,len(xDat)):
        A = (1.0/sqrt(2*math.pi*sigma**2))
        B = -((xDat[i]-mean)**2)/(2.0*sigma**2)
        fx = (hmax*(A*math.exp(B)))+hmin
        
        ci += ((yDat[i] - fx)**2)/ (yError[i]**2)
    return ci           

# Program for genorating a gaussian distribution         
def genGauss(hmax, hmin, mean, sigma, x):
    y = [] 
    for i in range(0,len(x)):
        A = (1.0/sqrt(2*math.pi*sigma**2.0))
        B = -((x[i]-mean)**2.0)/(2.0*sigma**2.0)
        temp = (hmax*A*math.exp(B))
        temp += hmin
        y.append(temp)
    return y

# for fitting a second order polynomial
def polyfitter(x,*args):    
    xDat = args[0]
    yDat = args[1]
    yError = 1.0
    A = x[0]
    B = x[1]
    C = x[2]
    ci = 0.0
    for i in range(0,len(xDat)):
        fx = A*xDat[i]**2 + B*xDat[i] + C
        ci += ((yDat[i] - fx)**2)/ (yError**2)
    return ci
 # for genorating a second order polynomial     
def polyGen(x,xval):
     xGen = range(0,int(max(xval)))
     y = []
     for i in range(0,len(xGen)):
         fx = []
         fx = fx = x[0]*xGen[i]**2 + x[1]*xGen[i] + x[2]
         y = append(y,fx)
     return xGen,y

#################   Source extraction ####################################
#This function takes the raw funtools output and then turn it into a flux
def funtoflux(filename, band):
    if int(band) == 250:
        x = 59.1
    elif int(band) == 350:
        x = 133.0
    elif int(band) == 500:
        x = 63.0
    elif int(band) == 160 or int(band) == 70 or int(band) == 100:
        x = 1000.0
    else:
        raise 'no band set'
    output = ["flux " + str(band) +"um" , "pixels"]
    output = array(output,ndmin=2)  
    output = output.reshape(1,2)
    funtoolsfile = open(pj(folder, filename), 'r')
    count = 0
    for line in funtoolsfile.readlines():
        count += 1
        try: 
            if line[3] == "1" and count > 16 and line[2] == " ":
                firstline = count -1
            else:
                continue
        except:
            continue
    count = 0    
    funtoolsfile.close()
    funtoolsfile = open(pj(folder, filename), 'r')    
    for row in funtoolsfile.readlines():
        count += 1
        info = row.split()
        try:        
            if count > firstline:
                selection = array((str(x*float(info[1])), str(info[2])))
                selection = selection.reshape(1,len(selection)) 
                output = append(output,selection, axis=0)
        except:
            continue 
    return output

#taking values for flux in jy or mjy and then taking away th ebackground and working out s/n
# this program requires an array to be defined globally call errorVal that has all the values 
# of sigma for combinded signel to noise.
def backgroundsubtraction(flux,backgrounds,band):
    table = array([(flux[0,0] + " backgroundsubtracted"),flux[0,0],"background", "detection threshold",("detected"+band)],ndmin=2)
    if band == "160":
        errorVal = errorValues[:,0]     
    elif band == "250":
        errorVal = errorValues[:,1]  
    elif band == "350":
        errorVal = errorValues[:,2]  
    elif band == "500":
        errorVal = errorValues[:,3]  
    else:
        raise 'no noise information present'
    count = 0  
    detected = 0
    for i in range(1,len(flux)):
        d =[]
        count += 1
        rawflux = float(flux[i,0])
        try:
            back = (float(backgrounds[i,0])/float(backgrounds[i,1]))*float(flux[i,1])
        except:
            back = 0.0
        fluxSub = rawflux - back
        sigmaFlux = errorVal[0]*float(flux[i,1])**2 + errorVal[1]*float(flux[i,1]) + errorVal[2]
        if fluxSub >= 5.0*sigmaFlux:
            detected += 1
            d = 1
        else:
            d = 0
        selection = array([fluxSub, rawflux, back, (5.0*sigmaFlux), d],ndmin=2)
        table = append(table,selection, axis=0)
    print band," detections = ",detected
    return table 

# for a quick way of finding the noise through working out 3 sigma of the stdev and then clipping 
# and then waiting untill it conveges.
def noisefinder(flux):
    a = flux[1:,0]
    b = array(a,dtype=float)
    sigma = std(b)
    newsigma = sigma*2.0
    count = 0
    while newsigma/sigma > 0.999:    
        sigma = std(b)
        meanF =[0.0]
        meanF = mean(b)
        selection = where((b < (meanF + 3*sigma))&(b > (meanF - 3*sigma)))
        b = b[selection]
        count += 1
        newsigma = std(b)
    return newsigma 
    

