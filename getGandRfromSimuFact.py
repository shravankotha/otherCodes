import sys

def main():
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
    print(listElementsWithInterface)
    

def findInterfacesSolidLiquidInAnElement(listNodesConnectedElement,listTemperaturesElement,listCoordinatesNodalElement):
    listXi,listEta,listZeta=coordsNaturalElementHex()
    nDivisions=4
    for iDivision in range(0,nDivisions):
        list
        for iCube in range(0,listCubes)
    

def 
    
    
def evaluateShapeFunctions(coordinatesNatural):
    listXi,listEta,listZeta=coordsNaturalElementHex()
    xi=coordinatesNatural[0]
    eta=coordinatesNatural[1]
    zeta=coordinatesNatural[2]
    for iNode in range(0,8):
        shapeFunctions[iNode]=(1/8)*(1+listXi[iNode]*xi)*(1+listEta[iNode]*eta)*(1+listZeta[iNode]*zeta)
        
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
        coordinatesNewElementsHex[0][iDimension].append([coordinatesNodesOld[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0]])
        coordinatesNewElementsHex[1][iDimension].append([coordinatesNodesNew[iDimension][0],coordinatesNodesOld[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0]])
        coordinatesNewElementsHex[2][iDimension].append([coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesOld[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0]])
        coordinatesNewElementsHex[3][iDimension].append([coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesOld[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0]])
        coordinatesNewElementsHex[4][iDimension].append([coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesOld[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0]])
        coordinatesNewElementsHex[5][iDimension].append([coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesOld[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0]])
        coordinatesNewElementsHex[6][iDimension].append([coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesOld[iDimension][0],coordinatesNodesNew[iDimension][0]])
        coordinatesNewElementsHex[7][iDimension].append([coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesNew[iDimension][0],coordinatesNodesOld[iDimension][0]])
    
    return coordinatesNodesNew, coordinatesNewElementsHex
    
def dShapeFunction_dNaturalCoords():
    listXi,listEta,listZeta=coordsNaturalElementHex()
    for iNode in range(0,8):        
        dN_dXi[iNode] = (1/8)*listXi[iNode]*(1+listEta[iNode]*eta)*(1+listZeta[iNode]*zeta)
        dN_dEta[iNode] = (1/8)*(1+listXi[iNode]*xi)*listEta[iNode]*(1+listZeta[iNode]*zeta)
        dN_dZeta[iNode] = (1/8)*(1+listXi[iNode]*xi)*(1+listEta[iNode]*eta)*listZeta[iNode]
               
def coordsNaturalElementHex():
    listXi=[]
    listEta=[]
    listZeta=[]
    
    listXi[0] = -1
    listXi[1] = +1
    listXi[2] = +1
    listXi[3] = -1
    listXi[4] = -1
    listXi[5] = +1
    listXi[6] = +1
    listXi[7] = -1
    
    listEta[0] = -1
    listEta[1] = -1
    listEta[2] = +1
    listEta[3] = +1
    listEta[4] = -1
    listEta[5] = -1
    listEta[6] = +1
    listEta[7] = +1    
    
    listZeta[0] = -1
    listZeta[1] = -1
    listZeta[2] = -1
    listZeta[3] = -1
    listZeta[4] = +1
    listZeta[5] = +1
    listZeta[6] = +1
    listZeta[7] = +1
    
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
            if listTemperatures[idNode] < temperatureLiquidus:
                isFoundNodeWithLowerThanLiquidus="true"
                idNodeWithLowerThanLiquidus=idNode
            elif listTemperatures[idNode] >= temperatureLiquidus:
                isFoundNodeWithHigherThanLiquidus="true"
                idNodeWithHigherThanLiquidus=idNode
            else:
                raise ValueError("Error: Something wrong with the temperature value")
                
            if isFoundNodeWithLowerThanLiquidus == "true" and isFoundNodeWithHigherThanLiquidus == "true":
                listElementsWithInterface.append(iElement)
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
