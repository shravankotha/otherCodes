import sys
import statistics as stats
import random as rand
import time
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

def main():
    nSamples=5
    nArguments=len(sys.argv)
    if nArguments != 4:
        print('Error: Three command line arguments are expected -- (1) unv file name including path (2) liquidus temperature in Kelvin (3) Temperature units (C or K)')
        return
    nameFile = sys.argv[1]
    temperatureLiquidus = float(sys.argv[2])
    unitsTemperature = sys.argv[3]
    
    
    listNodes,listCoordinatesX,listCoordinatesY,listCoordinatesZ,listConnectivity,listTemperatures = parseUnvFile(nameFile)
    listElementsWithInterface=findInterfaceSolidLiquid(listConnectivity, \
                                                       listTemperatures, \
                                                       unitsTemperature, \
                                                       temperatureLiquidus)
    
    listTemperaturesPointsOnInterface=[]
    listCoordinatesPointsOnInterface=[] 
    for iDimension in range(0,3):
        listCoordinatesPointsOnInterface.append([])
    timeStart=time.time()
    for iElement in range(0,500): #len(listElementsWithInterface)):
        print("Processing element: ",iElement," of ",str(len(listElementsWithInterface))," -- total time : ",str(time.time()-timeStart))
        listNodesConnectedElement= listConnectivity[listElementsWithInterface[iElement]-1]
        listTemperaturesElement = [listTemperatures[iNode-1] for iNode in listNodesConnectedElement]
        listCoordinatesNodalElement=[]
        for iDimension in range(0,3):
            listCoordinatesNodalElement.append([])  
            for iNode in range(0,len(listNodesConnectedElement)):
                if iDimension == 0:
                    listCoordinatesNodalElement[iDimension].append(listCoordinatesX[listNodesConnectedElement[iNode]-1])
                elif iDimension == 1:
                    listCoordinatesNodalElement[iDimension].append(listCoordinatesY[listNodesConnectedElement[iNode]-1])
                elif iDimension == 2:
                    listCoordinatesNodalElement[iDimension].append(listCoordinatesZ[listNodesConnectedElement[iNode]-1])
        
        coordinatesNodesNewWithInterface, temperatureNodesNewWithInterface = findInterfacesSolidLiquidInAnElement(listNodesConnectedElement,listTemperaturesElement,listCoordinatesNodalElement,temperatureLiquidus)
        
        #  evaluate the centroid of each cube, temperature at the centroid and gradient at the centroid
        #listTemperaturesPointsOnInterface.append([])
        #listCoordinatesPointsOnInterface.append([])
        listCubeIDs = rand.sample(range(len(coordinatesNodesNewWithInterface)),min([len(coordinatesNodesNewWithInterface),nSamples]))
        for iCube in range(0,len(listCubeIDs)):
            idCube = listCubeIDs[iCube]
            #listCoordinatesPointsOnInterface[iElement].append([])
            listTemperaturesCubeCorners=[]
            listCartesianCoordinates=[]
            for iDimension in range(0,3):
                listCartesianCoordinates.append([])
            for iNode in range(0,len(coordinatesNodesNewWithInterface[0][0])):
                coordinatesNatural=[]
                for iDimension in range(0,3):
                    coordinatesNatural.append(coordinatesNodesNewWithInterface[idCube][iDimension][iNode])
                listTemperaturesCubeCorners.append(getQuantityInsideTheElement(listTemperaturesElement,coordinatesNatural))
                for iDimension in range(0,3):
                    listCartesianCoordinates[iDimension].append(getQuantityInsideTheElement(listCoordinatesNodalElement[iDimension],coordinatesNatural))
        
            listTemperaturesPointsOnInterface.append(stats.mean(listTemperaturesCubeCorners))
            for iDimension in range(0,3):
                listCoordinatesPointsOnInterface[iDimension].append(stats.mean(listCartesianCoordinates[iDimension]))
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    #
    ## Data for a three-dimensional line
    #zline = np.linspace(0, 15, 1000)
    #xline = np.sin(zline)
    #yline = np.cos(zline)
    #ax.plot3D(xline, yline, zline, 'gray')
    #
    ## Data for three-dimensional scattered points
    #zdata = 15 * np.random.random(100)
    #xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
    #ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
    ax.scatter3D(listCoordinatesPointsOnInterface[0], listCoordinatesPointsOnInterface[1], listCoordinatesPointsOnInterface[2], color='red');
    plt.show()
    
    print('Finished')
    

    

def findInterfacesSolidLiquidInAnElement(listNodesConnectedElement,listTemperaturesElement,listCoordinatesNodalElement,temperatureLiquidus):
    listXi,listEta,listZeta=coordsNaturalElementHex()
    coordinatesNodesOld=[]
    coordinatesNodesOld.append([listXi,listEta,listZeta])
    nDivisions=6
    for iDivision in range(0,nDivisions):
        coordinatesNodesNewWithInterface = []
        temperatureNodesNewWithInterface = []
        for iCube in range(0,len(coordinatesNodesOld)):
            coordinatesNodesNew, coordinatesNewElementsHex = partitionCube(coordinatesNodesOld[iCube])
            for iCubeNew in range(0,len(coordinatesNewElementsHex)):
                coordinatesNodeElemental = coordinatesNewElementsHex[iCubeNew]
                isFoundNodeWithLowerThanLiquidus = "false"
                isFoundNodeWithHigherThanLiquidus = "false"
                listTemperatures=[]
                for iNode in range(0,len(coordinatesNodeElemental[0])):
                    coordinatesNatural = [coordinatesNodeElemental[0][iNode],coordinatesNodeElemental[1][iNode],coordinatesNodeElemental[2][iNode]]
                    listTemperatures.append(getQuantityInsideTheElement(listTemperaturesElement,coordinatesNatural))
                    if listTemperatures[iNode] < temperatureLiquidus:
                        isFoundNodeWithLowerThanLiquidus = "true"
                    elif listTemperatures[iNode] >= temperatureLiquidus:
                        isFoundNodeWithHigherThanLiquidus = "true"
                    
                if isFoundNodeWithLowerThanLiquidus == "true" and isFoundNodeWithHigherThanLiquidus == "true":
                    coordinatesNodesNewWithInterface.append(coordinatesNodeElemental)
                    temperatureNodesNewWithInterface.append(listTemperatures)
        coordinatesNodesOld = coordinatesNodesNewWithInterface      

    return coordinatesNodesNewWithInterface, temperatureNodesNewWithInterface
    
def getQuantityInsideTheElement(quantitiesNodal,coordinatesNatural):
    shapeFunctions=evaluateShapeFunctions(coordinatesNatural)
    quantity = 0
    for iNode in range(0,len(quantitiesNodal)):
        quantity = quantity + shapeFunctions[iNode]*quantitiesNodal[iNode]
        
    return quantity
    
def evaluateShapeFunctions(coordinatesNatural):
    listXi,listEta,listZeta=coordsNaturalElementHex()
    xi=coordinatesNatural[0]
    eta=coordinatesNatural[1]
    zeta=coordinatesNatural[2]
    shapeFunctions=[]
    for iNode in range(0,8):
        shapeFunctions.append((1/8)*(1+listXi[iNode]*xi)*(1+listEta[iNode]*eta)*(1+listZeta[iNode]*zeta))
        
    return shapeFunctions    
        
    
def partitionCube(coordinatesNodesOld):
    coordinatesNodesNew=[]    
    for iDimension in range(0,3):
        coordinatesNodesNew.append([])
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][0]+coordinatesNodesOld[iDimension][1]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][1]+coordinatesNodesOld[iDimension][2]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][2]+coordinatesNodesOld[iDimension][3]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][3]+coordinatesNodesOld[iDimension][0]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][4]+coordinatesNodesOld[iDimension][5]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][5]+coordinatesNodesOld[iDimension][6]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][6]+coordinatesNodesOld[iDimension][7]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][7]+coordinatesNodesOld[iDimension][4]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][1]+coordinatesNodesOld[iDimension][5]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][2]+coordinatesNodesOld[iDimension][6]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][3]+coordinatesNodesOld[iDimension][7]))
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][0]+coordinatesNodesOld[iDimension][4]))
        
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][0]+coordinatesNodesOld[iDimension][1]+coordinatesNodesOld[iDimension][2]+coordinatesNodesOld[iDimension][3]))
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][4]+coordinatesNodesOld[iDimension][5]+coordinatesNodesOld[iDimension][6]+coordinatesNodesOld[iDimension][7]))
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][0]+coordinatesNodesOld[iDimension][4]+coordinatesNodesOld[iDimension][5]+coordinatesNodesOld[iDimension][1]))
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][3]+coordinatesNodesOld[iDimension][7]+coordinatesNodesOld[iDimension][6]+coordinatesNodesOld[iDimension][2]))
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][0]+coordinatesNodesOld[iDimension][3]+coordinatesNodesOld[iDimension][7]+coordinatesNodesOld[iDimension][4]))
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][1]+coordinatesNodesOld[iDimension][2]+coordinatesNodesOld[iDimension][6]+coordinatesNodesOld[iDimension][5]))
        
        coordinatesNodesNew[iDimension].append((1/8)*(coordinatesNodesOld[iDimension][0]+coordinatesNodesOld[iDimension][1]+coordinatesNodesOld[iDimension][2]+coordinatesNodesOld[iDimension][3] + \
                                                  coordinatesNodesOld[iDimension][4]+coordinatesNodesOld[iDimension][5]+coordinatesNodesOld[iDimension][6]+coordinatesNodesOld[iDimension][7]))
    
    coordinatesNewElementsHex=[]
    for iCube in range(0,8):
        coordinatesNewElementsHex.append([])
        for iDimension in range(0,3):
            coordinatesNewElementsHex[iCube].append([])
        
    for iDimension in range(0,3):
        coordinatesNewElementsHex[0][iDimension].append(coordinatesNodesOld[iDimension][0])
        coordinatesNewElementsHex[0][iDimension].append(coordinatesNodesNew[iDimension][0])
        coordinatesNewElementsHex[0][iDimension].append(coordinatesNodesNew[iDimension][12])
        coordinatesNewElementsHex[0][iDimension].append(coordinatesNodesNew[iDimension][3])
        coordinatesNewElementsHex[0][iDimension].append(coordinatesNodesNew[iDimension][11])
        coordinatesNewElementsHex[0][iDimension].append(coordinatesNodesNew[iDimension][14])
        coordinatesNewElementsHex[0][iDimension].append(coordinatesNodesNew[iDimension][18])
        coordinatesNewElementsHex[0][iDimension].append(coordinatesNodesNew[iDimension][16])
        coordinatesNewElementsHex[1][iDimension].append(coordinatesNodesNew[iDimension][0])
        coordinatesNewElementsHex[1][iDimension].append(coordinatesNodesOld[iDimension][1])
        coordinatesNewElementsHex[1][iDimension].append(coordinatesNodesNew[iDimension][1])
        coordinatesNewElementsHex[1][iDimension].append(coordinatesNodesNew[iDimension][12])
        coordinatesNewElementsHex[1][iDimension].append(coordinatesNodesNew[iDimension][14])
        coordinatesNewElementsHex[1][iDimension].append(coordinatesNodesNew[iDimension][8])
        coordinatesNewElementsHex[1][iDimension].append(coordinatesNodesNew[iDimension][17])
        coordinatesNewElementsHex[1][iDimension].append(coordinatesNodesNew[iDimension][18])
        coordinatesNewElementsHex[2][iDimension].append(coordinatesNodesNew[iDimension][12])
        coordinatesNewElementsHex[2][iDimension].append(coordinatesNodesNew[iDimension][1])
        coordinatesNewElementsHex[2][iDimension].append(coordinatesNodesOld[iDimension][2])
        coordinatesNewElementsHex[2][iDimension].append(coordinatesNodesNew[iDimension][2])
        coordinatesNewElementsHex[2][iDimension].append(coordinatesNodesNew[iDimension][18])
        coordinatesNewElementsHex[2][iDimension].append(coordinatesNodesNew[iDimension][17])
        coordinatesNewElementsHex[2][iDimension].append(coordinatesNodesNew[iDimension][9])
        coordinatesNewElementsHex[2][iDimension].append(coordinatesNodesNew[iDimension][15])
        coordinatesNewElementsHex[3][iDimension].append(coordinatesNodesNew[iDimension][3])
        coordinatesNewElementsHex[3][iDimension].append(coordinatesNodesNew[iDimension][3])
        coordinatesNewElementsHex[3][iDimension].append(coordinatesNodesNew[iDimension][2])
        coordinatesNewElementsHex[3][iDimension].append(coordinatesNodesOld[iDimension][3])
        coordinatesNewElementsHex[3][iDimension].append(coordinatesNodesNew[iDimension][16])
        coordinatesNewElementsHex[3][iDimension].append(coordinatesNodesNew[iDimension][18])
        coordinatesNewElementsHex[3][iDimension].append(coordinatesNodesNew[iDimension][15])
        coordinatesNewElementsHex[3][iDimension].append(coordinatesNodesNew[iDimension][10])
        coordinatesNewElementsHex[4][iDimension].append(coordinatesNodesNew[iDimension][11])
        coordinatesNewElementsHex[4][iDimension].append(coordinatesNodesNew[iDimension][14])
        coordinatesNewElementsHex[4][iDimension].append(coordinatesNodesNew[iDimension][18])
        coordinatesNewElementsHex[4][iDimension].append(coordinatesNodesNew[iDimension][16])
        coordinatesNewElementsHex[4][iDimension].append(coordinatesNodesOld[iDimension][4])
        coordinatesNewElementsHex[4][iDimension].append(coordinatesNodesNew[iDimension][4])
        coordinatesNewElementsHex[4][iDimension].append(coordinatesNodesNew[iDimension][13])
        coordinatesNewElementsHex[4][iDimension].append(coordinatesNodesNew[iDimension][7])
        coordinatesNewElementsHex[5][iDimension].append(coordinatesNodesNew[iDimension][14])
        coordinatesNewElementsHex[5][iDimension].append(coordinatesNodesNew[iDimension][8])
        coordinatesNewElementsHex[5][iDimension].append(coordinatesNodesNew[iDimension][17])
        coordinatesNewElementsHex[5][iDimension].append(coordinatesNodesNew[iDimension][18])
        coordinatesNewElementsHex[5][iDimension].append(coordinatesNodesNew[iDimension][4])
        coordinatesNewElementsHex[5][iDimension].append(coordinatesNodesOld[iDimension][5])
        coordinatesNewElementsHex[5][iDimension].append(coordinatesNodesNew[iDimension][5])
        coordinatesNewElementsHex[5][iDimension].append(coordinatesNodesNew[iDimension][13])
        coordinatesNewElementsHex[6][iDimension].append(coordinatesNodesNew[iDimension][18])
        coordinatesNewElementsHex[6][iDimension].append(coordinatesNodesNew[iDimension][17])
        coordinatesNewElementsHex[6][iDimension].append(coordinatesNodesNew[iDimension][9])
        coordinatesNewElementsHex[6][iDimension].append(coordinatesNodesNew[iDimension][15])
        coordinatesNewElementsHex[6][iDimension].append(coordinatesNodesNew[iDimension][13])
        coordinatesNewElementsHex[6][iDimension].append(coordinatesNodesNew[iDimension][5])
        coordinatesNewElementsHex[6][iDimension].append(coordinatesNodesOld[iDimension][6])
        coordinatesNewElementsHex[6][iDimension].append(coordinatesNodesNew[iDimension][6])
        coordinatesNewElementsHex[7][iDimension].append(coordinatesNodesNew[iDimension][16])
        coordinatesNewElementsHex[7][iDimension].append(coordinatesNodesNew[iDimension][18])
        coordinatesNewElementsHex[7][iDimension].append(coordinatesNodesNew[iDimension][15])
        coordinatesNewElementsHex[7][iDimension].append(coordinatesNodesNew[iDimension][10])
        coordinatesNewElementsHex[7][iDimension].append(coordinatesNodesNew[iDimension][7])
        coordinatesNewElementsHex[7][iDimension].append(coordinatesNodesNew[iDimension][13])
        coordinatesNewElementsHex[7][iDimension].append(coordinatesNodesNew[iDimension][6])
        coordinatesNewElementsHex[7][iDimension].append(coordinatesNodesOld[iDimension][7])
    
    return coordinatesNodesNew, coordinatesNewElementsHex
    
def dShapeFunction_dNaturalCoords():
    listXi,listEta,listZeta=coordsNaturalElementHex()
    dN_dXi=[]
    dN_dEta=[]
    dN_dZeta=[]
    for iNode in range(0,8):        
        dN_dXi.append((1/8)*listXi[iNode]*(1+listEta[iNode]*eta)*(1+listZeta[iNode]*zeta))
        dN_dEta.append((1/8)*(1+listXi[iNode]*xi)*listEta[iNode]*(1+listZeta[iNode]*zeta))
        dN_dZeta.append((1/8)*(1+listXi[iNode]*xi)*(1+listEta[iNode]*eta)*listZeta[iNode])
               
def coordsNaturalElementHex():
    listXi=[]
    listEta=[]
    listZeta=[]
    
    listXi.append(-1)
    listXi.append(+1)
    listXi.append(+1)
    listXi.append(-1)
    listXi.append(-1)
    listXi.append(+1)
    listXi.append(+1)
    listXi.append(-1)
    
    listEta.append(-1)
    listEta.append(-1)
    listEta.append(+1)
    listEta.append(+1)
    listEta.append(-1)
    listEta.append(-1)
    listEta.append(+1)
    listEta.append(+1)    
    
    listZeta.append(-1)
    listZeta.append(-1)
    listZeta.append(-1)
    listZeta.append(-1)
    listZeta.append(+1)
    listZeta.append(+1)
    listZeta.append(+1)
    listZeta.append(+1)
    
    return listXi, listEta, listZeta



def findInterfaceSolidLiquid(listConnectivity, \
                             listTemperatures, \
                             unitsTemperature, \
                             temperatureLiquidus):
    listElementsWithInterface=[]
    for iElement in range(0,len(listConnectivity)):
        listConnectedNodes=listConnectivity[iElement]
        isFoundNodeWithLowerThanLiquidus="false"
        isFoundNodeWithHigherThanLiquidus="false"
        for iNode in range(0,len(listConnectedNodes)):
            idNode=listConnectedNodes[iNode]-1
            temperature=listTemperatures[idNode]
            if temperature < temperatureLiquidus:
                isFoundNodeWithLowerThanLiquidus="true"
                idNodeWithLowerThanLiquidus=idNode
            elif temperature >= temperatureLiquidus:
                isFoundNodeWithHigherThanLiquidus="true"
                idNodeWithHigherThanLiquidus=idNode
            else:
                raise ValueError("Error: Something wrong with the temperature value")
                
            if isFoundNodeWithLowerThanLiquidus == "true" and isFoundNodeWithHigherThanLiquidus == "true":
                listElementsWithInterface.append(iElement+1)
                break
                
    return listElementsWithInterface
                
   
            
            
def multiplyMatrices(A,B):
    if len(A[0]) != len(B):
        raise ValuError("Error: The dimensions of the matrices are inconsistent for multiplication")
    C=[]
    for ii in range(0,len(A)):
        C.append([])
        C[ii] = [0]*len(B[0])
        for jj in range(0,len(B[0])):            
            for kk in range(0,len(A[ii])):
                C[ii][jj] = C[ii][jj] + A[ii][kk]*B[kk][jj]
    return C            
        
            
    
def parseUnvFile(nameFile):
    listNodesTemperature=[]
    listNodesNodalCoordinates=[]
    listIdElementsConnectivity=[] 
    with open(nameFile) as fileCurrent:
        dataTemp=fileCurrent.readlines()
        isFoundKeywordNodes='false'
        isFoundKeywordConnectivity='false'
        isFoundKeywordTemp='false'
        shouldReadNodes='false'
        shouldKeepReadingNodes='true'
        shouldKeepReadingConnectivity='true'
        count=0
        for iData in range(0,len(dataTemp)):
            # Reading the nodal coordinates
            if str(dataTemp[iData].strip())=="Contained bodies:":
                isFoundKeywordNodes='true'
                count=0 
                
            if isFoundKeywordNodes=='true' and shouldKeepReadingNodes=='true':                
                count=count+1            
                if count>=6:
                    listNodesNodalCoordinates.append(dataTemp[iData].strip().split())
                if str(dataTemp[iData].strip())=="-1" and str(dataTemp[iData+1].strip())=="-1":
                    shouldKeepReadingNodes='false'
                
            # Reading element connectivity
            if str(dataTemp[iData].strip())=="2412" and str(dataTemp[iData-1].strip())=="-1":
                isFoundKeywordConnectivity='true'
                count=0
            if isFoundKeywordConnectivity=='true' and shouldKeepReadingConnectivity=='true':
                count=count+1
                if count>1:
                    listIdElementsConnectivity.append(dataTemp[iData].strip().split())
                if str(dataTemp[iData].strip())=="-1" and str(dataTemp[iData+1].strip())=="-1":
                    shouldKeepReadingConnectivity='false'                    
                    
            # Reading the temperature data
            if str(dataTemp[iData].strip())=="Temperature":
                isFoundKeywordTemp='true'
                count=0     
                
            if isFoundKeywordTemp=='true':
                if count <= 9:    # skipping 11 lines after Temperature keyword is found
                    count=count+1
                    continue
                else:
                    if iData!=len(dataTemp)-1:
                        listNodesTemperature.append(dataTemp[iData].strip())        
                        
        # Do some checks and write coordinates and temperatures into lists
        if isFoundKeywordNodes=='false' or isFoundKeywordConnectivity=='false' or isFoundKeywordTemp=='false':
            sys.exit("Error: Nodal coords or connectivity or Temperature keyword/information is not found in the file")
        else:    
            listNodes=[]
            listCoordinatesX=[]
            listCoordinatesY=[]
            listCoordinatesZ=[]
            listElements=[]
            listConnectivity=[]
            # Place nodal coordinates into arrays
            for i in range(0,len(listNodesNodalCoordinates)-1,2):
                listNodes.append(int(listNodesNodalCoordinates[i][0]))
            for i in range(1,len(listNodesNodalCoordinates)-1,2):
                listCoordinatesX.append(float(listNodesNodalCoordinates[i][0]))
                listCoordinatesY.append(float(listNodesNodalCoordinates[i][1]))
                listCoordinatesZ.append(float(listNodesNodalCoordinates[i][2]))
            # Place element connectivity into arrays
            for i in range(0,len(listIdElementsConnectivity)-1,2):
                listElements.append(int(listIdElementsConnectivity[i][0]))
            count=0    
            for i in range(1,len(listIdElementsConnectivity)-1,2):
                count=count+1
                listConnectivity.append([])
                for j in range(0,8):                    
                    listConnectivity[count-1].append(int(listIdElementsConnectivity[i][j]))
            # Place temperature into arrays
            listNodesTempData=[int(i) for i in listNodesTemperature[0:len(listNodesTemperature):2]]
            listTemperatures=[float(i) for i in listNodesTemperature[1:len(listNodesTemperature):2]]

    return listNodes,listCoordinatesX,listCoordinatesY,listCoordinatesZ,listConnectivity,listTemperatures
    
if __name__ == "__main__":
    main()
