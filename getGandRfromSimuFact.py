import sys
import os
import statistics as stats
import random as rand
import time
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
import numpy as np
import math
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from interpolate2D import interpolate2D
from partitionCube import partitionCube
from hexElement import dShapeFunction_dCaterianCoords, getJacobianHexElement
from hexElement import dShapeFunction_dNaturalCoords, coordsNaturalElementHex
from hexElement import evaluateShapeFunctions, getQuantityInsideTheElement
from multiplyMatrices import multiplyMatrices
from parseUnvFile_SimuFact import parseUnvFile
from meltpoolCalculationsFromFE import findElementsContainingInterfaceSolidLiquid, findInterfacesSolidLiquidInAnElement
from plotMeltpool4D import plotMeltpoolScatter4D

def main():

    methodInterpolation = 'linear'
    
    nPoints = 50
    
    offsetTol = 0.005
    
    tolerance = 1E-6   
    
    timeStart = time.time()
    
    nArguments = len(sys.argv)
    
    if nArguments != 7:
    
        print('Error: Five command line arguments are expected -- (1) unv file name including path      \
                                                                  (2) liquidus temperature in Kelvin    \
                                                                  (3) Temperature units (C or K)        \
                                                                  (4) Laser direction (X or Y or Z)     \
                                                                  (5) Laser speed                       \
                                                                  (6) Melt pool normal direction')
                                                                  
        return
        
    nameFile = sys.argv[1]
    
    temperatureLiquidus = float(sys.argv[2])
    
    unitsTemperature = sys.argv[3]
    
    directionScanning = sys.argv[4]
    
    speedLaser = float(sys.argv[5])
    
    directionMeltPoolNormal = sys.argv[6]
    
    iDimensionScanning = {"X":0, "Y":1, "Z":2}
    
    # -------------------------------------- Parse the unv file and obtain mesh information
    
    listNodes, listCoordinatesX, listCoordinatesY, listCoordinatesZ, listConnectivity, listTemperatures = parseUnvFile(nameFile)
    
    listElementIDs = [ii+1 for ii in range(0,len(listConnectivity))]
    
    # -------------------------------------- Find the solid-liquid interface
    
    listElementsWithInterface = findElementsContainingInterfaceSolidLiquid(listElementIDs,
                                                                           listConnectivity,
                                                                           listTemperatures,
                                                                           unitsTemperature,
                                                                           temperatureLiquidus)
    
    listTemperaturesPointsOnInterface, listCoordinatesCartesianPointsOnInterface, listCoordinatesNaturalPointsOnInterface, \
            listVectorGradient, listMagnitudeGradient, listDendriteGrowthRate, listGrowthRate, listCoolingRate = ([] for ii in range(8))
    
    for iDimension in range(0, 3):
    
        listCoordinatesCartesianPointsOnInterface.append([])
        
        listCoordinatesNaturalPointsOnInterface.append([])
        
        listVectorGradient.append([])
    
    # -------------------------------------- Loop over all the interface elements and compute the interface points, thermal gradient, growth rate and cooling rate
    
    for iElement in range(0, len(listElementsWithInterface), 1):
    
        print("Processing element: ", iElement, " of ", str(len(listElementsWithInterface)), " -- elapsed time : ", str(time.time()-timeStart))
        
        listNodesConnectedElement = listConnectivity[listElementsWithInterface[iElement] - 1]
        
        listTemperaturesElement = [listTemperatures[iNode - 1] for iNode in listNodesConnectedElement]
        
        listCoordinatesNodalElement = []
        
        for iDimension in range(0, 3):
        
            listCoordinatesNodalElement.append([])  
            
            for iNode in range(0, len(listNodesConnectedElement)):
            
                if iDimension == 0:
                
                    listCoordinatesNodalElement[iDimension].append(listCoordinatesX[listNodesConnectedElement[iNode] - 1])
                    
                elif iDimension == 1:
                
                    listCoordinatesNodalElement[iDimension].append(listCoordinatesY[listNodesConnectedElement[iNode] - 1])
                    
                elif iDimension == 2:
                
                    listCoordinatesNodalElement[iDimension].append(listCoordinatesZ[listNodesConnectedElement[iNode] - 1])
        
        coordinatesNodesNewWithInterface, temperatureNodesNewWithInterface = findInterfacesSolidLiquidInAnElement(listNodesConnectedElement,
                                                                                                                  listTemperaturesElement,
                                                                                                                  listCoordinatesNodalElement,
                                                                                                                  temperatureLiquidus)
        
        # -------------------------------------- evaluate the centroid of each cube, temperature at the centroid and gradient at the centroid
        
        for iCube in range(0, len(coordinatesNodesNewWithInterface)):
        
            listTemperaturesCubeCorners = []
            
            listCartesianCoordinates = []
            
            for iDimension in range(0, 3):
            
                listCartesianCoordinates.append([])
                
            for iNode in range(0, len(coordinatesNodesNewWithInterface[0][0])):
            
                coordinatesNatural = []
                
                for iDimension in range(0, 3):
                
                    coordinatesNatural.append(coordinatesNodesNewWithInterface[iCube][iDimension][iNode])
                    
                listTemperaturesCubeCorners.append(getQuantityInsideTheElement(listTemperaturesElement, coordinatesNatural))
                
                for iDimension in range(0, 3):
                
                    listCartesianCoordinates[iDimension].append(getQuantityInsideTheElement(listCoordinatesNodalElement[iDimension], coordinatesNatural))
        
            listTemperaturesPointsOnInterface.append(stats.mean(listTemperaturesCubeCorners))
            
            coordinatesNaturalCentroid = []
            
            for iDimension in range(0, 3):
            
                listCoordinatesCartesianPointsOnInterface[iDimension].append(stats.mean(listCartesianCoordinates[iDimension]))
                
                coordinatesNaturalCentroid.append(stats.mean(coordinatesNodesNewWithInterface[iCube][iDimension]))
                
                listCoordinatesNaturalPointsOnInterface[iDimension].append(coordinatesNaturalCentroid[iDimension])
                
            # -------------------------------------- Compute the gradient, growth rate
            
            dN_dxyz = dShapeFunction_dCaterianCoords(listCoordinatesNodalElement, coordinatesNaturalCentroid)
            
            magnitude = 0
            
            for iDimension in range(0, 3):
            
                gradient_ = 0
                
                for iNode in range(0, 8):
                
                    gradient_ = gradient_ + dN_dxyz[iNode][iDimension]*listTemperaturesElement[iNode]
                    
                magnitude = magnitude + gradient_**2
                
                listVectorGradient[iDimension].append(gradient_)
                
            listMagnitudeGradient.append(math.sqrt(magnitude))
            
            listDendriteGrowthRate.append(speedLaser*abs(listVectorGradient[iDimensionScanning[directionScanning]][len(listMagnitudeGradient) - 1])/max([abs(listVectorGradient[0][len(listMagnitudeGradient) - 1]),abs(listVectorGradient[1][len(listMagnitudeGradient) - 1]),abs(listVectorGradient[2][len(listMagnitudeGradient) - 1])]))   # from: https://www.sciencedirect.com/science/article/abs/pii/S1005030216300615
            
            listGrowthRate.append(speedLaser*abs(listVectorGradient[iDimensionScanning[directionScanning]][len(listMagnitudeGradient) - 1])/math.sqrt(magnitude))
            
            listCoolingRate.append(listMagnitudeGradient[len(listMagnitudeGradient) - 1]*listGrowthRate[len(listMagnitudeGradient) - 1])
    
    # -------------------------------------- sort the data points based on thermal gradient and write them to a file
    
    listIndices = np.argsort(np.array(listMagnitudeGradient))
    
    pathDir = os.path.splitdrive(nameFile)[0] + os.path.split(os.path.splitdrive(nameFile)[1])[0] + "/"
    
    out_path = pathDir + 'coords_Temp_G_R_Vd_coolingRate.out'
    
    with open(out_path, 'w') as file_out:
    
        file_out.write("coordsX     coordsY     coordsZ     Temp        thermalGradX        thermalGradY        thermalGradZ        thermalGradMagnitude        growthRate   dendriteGrowthRate   coolingRate \n")
        
        for ii in range(0,len(listIndices)):
        
            file_out.write("{0:25.10f}{1:25.10f}{2:25.10f}{3:25.10f}{4:25.10f}{5:25.10f}{6:25.10f}{7:25.10f}{8:25.10f}{9:25.10f}{10:25.10f}\n".format(listCoordinatesCartesianPointsOnInterface[0][listIndices[ii]], 
                                                                                         listCoordinatesCartesianPointsOnInterface[1][listIndices[ii]],
                                                                                         listCoordinatesCartesianPointsOnInterface[2][listIndices[ii]],
                                                                                         listTemperaturesPointsOnInterface[listIndices[ii]],
                                                                                         listVectorGradient[0][listIndices[ii]],
                                                                                         listVectorGradient[1][listIndices[ii]],
                                                                                         listVectorGradient[2][listIndices[ii]],
                                                                                         listMagnitudeGradient[listIndices[ii]],
                                                                                         listGrowthRate[listIndices[ii]],
                                                                                         listDendriteGrowthRate[listIndices[ii]],
                                                                                         listCoolingRate[listIndices[ii]]
                                                                                        )
                          )
                          
    file_out.close()
    
    # -------------------------------------- Find G and R for few important locations in the melt pool and write them to a file 
    
    listIndices = np.argsort(np.array(listMagnitudeGradient))
    
    lowestG = listMagnitudeGradient[listIndices[0]]
    
    moderateG = listMagnitudeGradient[listIndices[int(len(listIndices)/2)]]
    
    highestG = listMagnitudeGradient[listIndices[len(listIndices) - 1]]
    
    RforLowestG = listGrowthRate[listIndices[0]]
    
    RforModerateG = listGrowthRate[listIndices[int(len(listIndices)/2)]]
    
    RforHighestG = listGrowthRate[listIndices[len(listIndices) - 1]]
    
    VdforLowestG = listDendriteGrowthRate[listIndices[0]]
    
    VdforModerateG = listDendriteGrowthRate[listIndices[int(len(listIndices)/2)]]
    
    VdforHighestG = listDendriteGrowthRate[listIndices[len(listIndices) - 1]]
    
    
    # -------------------------------------- 
    listIndices = np.argsort(np.array(listGrowthRate))
    
    lowestR = listGrowthRate[listIndices[0]]
    
    moderateR = listGrowthRate[listIndices[int(len(listIndices)/2)]]
    
    highestR = listGrowthRate[listIndices[len(listIndices) - 1]]
    
    GforLowestR = listMagnitudeGradient[listIndices[0]]
    
    GforModerateR = listMagnitudeGradient[listIndices[int(len(listIndices)/2)]]
    
    GforHighestR = listMagnitudeGradient[listIndices[len(listIndices) - 1]]
    
    out_path = pathDir + 'G_R_forPFsimulations.out'
    
    with open(out_path, 'w') as file_out:
    
        file_out.write("Lowest G and corresponding R   : ")
        
        file_out.write("{0:25.10f}{1:25.10f}\n".format(lowestG, RforLowestG))
        
        file_out.write("Moderate G and corresponding R : ")
        
        file_out.write("{0:25.10f}{1:25.10f}\n".format(moderateG, RforModerateG))
        
        file_out.write("Highest G and corresponding R  : ")
        
        file_out.write("{0:25.10f}{1:25.10f}\n".format(highestG, RforHighestG))

        file_out.write("G corresponding to lowest R    : ")
        
        file_out.write("{0:25.10f}{1:25.10f}\n".format(GforLowestR, lowestR))
        
        file_out.write("G corresponding to moderate R  : ")
        
        file_out.write("{0:25.10f}{1:25.10f}\n".format(GforModerateR, moderateR))
        
        file_out.write("G corresponding to highest R   : ")
        
        file_out.write("{0:25.10f}{1:25.10f}\n".format(GforHighestR, highestR))
        
    file_out.close()
    
    # -------------------------------------- extract data along centerline of the meltpool
    
    coordinateMaxAlongMeltpool = max(listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]])
    
    coordinateMinAlongMeltpool = min(listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]])
    
    # find the other in-plane coordinate for the maximum and minimum above
    
    if directionScanning == "X" and directionMeltPoolNormal == "Y" or \
       directionScanning == "Y" and directionMeltPoolNormal == "X":
    
        directionTransverse = "Z"
        
    elif directionScanning == "X" and directionMeltPoolNormal == "Z" or \
       directionScanning == "Z" and directionMeltPoolNormal == "X":
    
        directionTransverse = "Y"

    elif directionScanning == "Y" and directionMeltPoolNormal == "Z" or \
       directionScanning == "Z" and directionMeltPoolNormal == "Y":
    
        directionTransverse = "X"        
    
    for iPoint in range(0,len(listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]])):
    
        if abs(coordinateMaxAlongMeltpool - listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]][iPoint]) <= tolerance*abs(coordinateMaxAlongMeltpool):
        
            idPoint = iPoint
            
            break
    
    coordinateTransverseForMaxCoordAlongMeltpool = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionTransverse]][idPoint]
    
    coordinateNormalForMaxCoordAlongMeltpool = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionMeltPoolNormal]][idPoint]
    
    # -------------------------------------- interpolate the melt pool normal coordinate, G, R and cooling rate along the centerline
        
    xMinToBeInterpolated = coordinateMinAlongMeltpool + offsetTol*(coordinateMaxAlongMeltpool - coordinateMinAlongMeltpool)
    
    xMaxToBeInterpolated = coordinateMaxAlongMeltpool - offsetTol*(coordinateMaxAlongMeltpool - coordinateMinAlongMeltpool)
    
    # interpolate melt pool depth
    coordsGridScanningDirection, coordsGridTransverseDirection, coordsGridNormalDirection = interpolate2D(                                          \
                                           listPointsXref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]],       \
                                           listPointsYref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionTransverse]],     \
                                           listPointsZref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionMeltPoolNormal]], \
                                           xMinToBeInterpolated = xMinToBeInterpolated,                                                             \
                                           xMaxToBeInterpolated = xMaxToBeInterpolated,                                                             \
                                           yMinToBeInterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                     \
                                           yMaxToBeinterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                     \
                                           nPoints = nPoints,                                                                                       \
                                           methodInterpolation = methodInterpolation)
                                            
    plt.plot(coordsGridScanningDirection[0], coordsGridNormalDirection[0])
    plt.title("Melt pool center line")
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    
    xlim_min = min(coordsGridScanningDirection[0]) - 6*offsetTol*(max(coordsGridScanningDirection[0])-min(coordsGridScanningDirection[0]))
    xlim_max = max(coordsGridScanningDirection[0]) + 6*offsetTol*(max(coordsGridScanningDirection[0])-min(coordsGridScanningDirection[0]))

    ylim_min = min(coordsGridNormalDirection[0]) - 6*offsetTol*(max(coordsGridNormalDirection[0])-min(coordsGridNormalDirection[0]))
    ylim_max = max(coordsGridNormalDirection[0]) + 6*offsetTol*(max(coordsGridNormalDirection[0])-min(coordsGridNormalDirection[0]))

    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    
    plt.show()
 
    # interpolate temperatures
    coordsGridScanningDirection, coordsGridTransverseDirection, gridTemperatures = interpolate2D(                                               \
                                           listPointsXref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]],   \
                                           listPointsYref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionTransverse]], \
                                           listPointsZref = listTemperaturesPointsOnInterface,                                                  \
                                           xMinToBeInterpolated = xMinToBeInterpolated,                                                         \
                                           xMaxToBeInterpolated = xMaxToBeInterpolated,                                                         \
                                           yMinToBeInterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           yMaxToBeinterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           nPoints = nPoints,                                                                                   \
                                           methodInterpolation = methodInterpolation)
                                            
    plt.figure()
    sc = plt.scatter(coordsGridScanningDirection[0], coordsGridNormalDirection[0], c = gridTemperatures[0], cmap = "jet")
    plt.colorbar()
    plt.title("Temperature (K) ")
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    plt.show()
 
    # interpolate thermal gradient
    coordsGridScanningDirection, coordsGridTransverseDirection, gridThermalGradient = interpolate2D(                                            \
                                           listPointsXref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]],   \
                                           listPointsYref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionTransverse]], \
                                           listPointsZref = listMagnitudeGradient,                                                              \
                                           xMinToBeInterpolated = xMinToBeInterpolated,                                                         \
                                           xMaxToBeInterpolated = xMaxToBeInterpolated,                                                         \
                                           yMinToBeInterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           yMaxToBeinterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           nPoints = nPoints,                                                                                   \
                                           methodInterpolation = methodInterpolation)
                                            
    plt.figure()
    sc = plt.scatter(coordsGridScanningDirection[0], coordsGridNormalDirection[0], c = gridThermalGradient[0], cmap = "jet")
    plt.colorbar()
    plt.title("Thermal gradient")
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    plt.show()
  
    # interpolate growth rate
    coordsGridScanningDirection, coordsGridTransverseDirection, gridGrowthRate = interpolate2D(                                                 \
                                           listPointsXref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]],   \
                                           listPointsYref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionTransverse]], \
                                           listPointsZref = listGrowthRate,                                                                     \
                                           xMinToBeInterpolated = xMinToBeInterpolated,                                                         \
                                           xMaxToBeInterpolated = xMaxToBeInterpolated,                                                         \
                                           yMinToBeInterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           yMaxToBeinterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           nPoints = nPoints,                                                                                   \
                                           methodInterpolation = methodInterpolation)
                                            
    plt.figure()
    sc = plt.scatter(coordsGridScanningDirection[0], coordsGridNormalDirection[0], c = gridGrowthRate[0], cmap = "jet")
    plt.colorbar()
    plt.title("Growth rate")
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    plt.show()    

    # interpolate dendrite growth rate
    coordsGridScanningDirection, coordsGridTransverseDirection, gridDendriteGrowthRate = interpolate2D(                                                 \
                                           listPointsXref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]],   \
                                           listPointsYref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionTransverse]], \
                                           listPointsZref = listDendriteGrowthRate,                                                             \
                                           xMinToBeInterpolated = xMinToBeInterpolated,                                                         \
                                           xMaxToBeInterpolated = xMaxToBeInterpolated,                                                         \
                                           yMinToBeInterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           yMaxToBeinterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           nPoints = nPoints,                                                                                   \
                                           methodInterpolation = methodInterpolation)
                                            
    plt.figure()
    sc = plt.scatter(coordsGridScanningDirection[0], coordsGridNormalDirection[0], c = gridDendriteGrowthRate[0], cmap = "jet")
    plt.colorbar()
    plt.title("Dendrite growth rate")
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    plt.show()
    
    # interpolate cooling rate
    coordsGridScanningDirection, coordsGridTransverseDirection, gridCoolingRate = interpolate2D(                                                \
                                           listPointsXref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionScanning]],   \
                                           listPointsYref = listCoordinatesCartesianPointsOnInterface[iDimensionScanning[directionTransverse]], \
                                           listPointsZref = listCoolingRate,                                                                    \
                                           xMinToBeInterpolated = xMinToBeInterpolated,                                                         \
                                           xMaxToBeInterpolated = xMaxToBeInterpolated,                                                         \
                                           yMinToBeInterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           yMaxToBeinterpolated = coordinateTransverseForMaxCoordAlongMeltpool,                                 \
                                           nPoints = nPoints,                                                                                   \
                                           methodInterpolation = methodInterpolation)
                                            
    plt.figure()
    sc = plt.scatter(coordsGridScanningDirection[0], coordsGridNormalDirection[0], c = gridCoolingRate[0], cmap = "jet")
    plt.colorbar()
    plt.title("Cooling rate")
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    plt.show()

    # write the data along the center line to file
    
    out_path = pathDir + 'coords_depth_Temp_G_R_Vd_coolingRate_alongCenterline.out'
    
    with open(out_path, 'w') as file_out:
    
        file_out.write("coordsLongitudinal     coordsTransverse     coordsNormal        depth          Temp        thermalGradMagnitude        growthRate   dendriteGrowthRate   coolingRate \n")
        
        for ii in range(0,len(coordsGridScanningDirection[0])):
        
            file_out.write("{0:25.10f}{1:25.10f}{2:25.10f}{3:25.10f}{4:25.10f}{5:25.10f}{6:25.10f}{7:25.10f}{8:25.10f}\n".format(coordsGridScanningDirection[0][ii], 
                                                                                         coordsGridTransverseDirection[0][ii],
                                                                                         coordsGridNormalDirection[0][ii],
                                                                                         abs(coordinateNormalForMaxCoordAlongMeltpool - coordsGridNormalDirection[0][ii]),
                                                                                         gridTemperatures[0][ii],                                                                                         
                                                                                         gridThermalGradient[0][ii],
                                                                                         gridGrowthRate[0][ii],
                                                                                         gridDendriteGrowthRate[0][ii],
                                                                                         gridCoolingRate[0][ii]
                                                                                        )
                          )
                          
    file_out.close()
    
    # -------------------------------------- 3D Plots 
    
    plotMeltpoolScatter4D(listCoordinatesCartesianPointsOnInterface[0],
                          listCoordinatesCartesianPointsOnInterface[1],
                          listCoordinatesCartesianPointsOnInterface[2],
                          listCoordinatesCartesianPointsOnInterface[2],
                          directionMeltPoolNormal,
                          "Melt Pool",
                          "X (m)",
                          "Y (m)",
                          "Z (m)",
                          shouldTurnOffGrid = False,
                          shouldTurnOffAxes = False)
                          
    plotMeltpoolScatter4D(listCoordinatesCartesianPointsOnInterface[0],
                          listCoordinatesCartesianPointsOnInterface[1],
                          listCoordinatesCartesianPointsOnInterface[2],
                          listMagnitudeGradient,
                          directionMeltPoolNormal,
                          "Thermal Gradient (G)",
                          "X (m)",
                          "Y (m)",
                          "Z (m)",
                          shouldTurnOffGrid = True,
                          shouldTurnOffAxes = True)
                          
    plotMeltpoolScatter4D(listCoordinatesCartesianPointsOnInterface[0],
                          listCoordinatesCartesianPointsOnInterface[1],
                          listCoordinatesCartesianPointsOnInterface[2],
                          listGrowthRate,
                          directionMeltPoolNormal,
                          "Growth Rate (R)",
                          "X (m)",
                          "Y (m)",
                          "Z (m)",
                          shouldTurnOffGrid = True,
                          shouldTurnOffAxes = True)

    plotMeltpoolScatter4D(listCoordinatesCartesianPointsOnInterface[0],
                          listCoordinatesCartesianPointsOnInterface[1],
                          listCoordinatesCartesianPointsOnInterface[2],
                          listDendriteGrowthRate,
                          directionMeltPoolNormal,
                          "Dendrite Growth Rate (R)",
                          "X (m)",
                          "Y (m)",
                          "Z (m)",
                          shouldTurnOffGrid = True,
                          shouldTurnOffAxes = True)
                          
    plotMeltpoolScatter4D(listCoordinatesCartesianPointsOnInterface[0],
                          listCoordinatesCartesianPointsOnInterface[1],
                          listCoordinatesCartesianPointsOnInterface[2],
                          listCoolingRate,
                          directionMeltPoolNormal,
                          "Cooling Rate",
                          "X (m)",
                          "Y (m)",
                          "Z (m)",
                          shouldTurnOffGrid = True,
                          shouldTurnOffAxes = True)

    # scatter plot the G vs R
    #plotScatter2D(listMagnitudeGradient,listGrowthRate,"G versus R","G (K/m)","R (m/s)")
    #plotMeltpoolSurface4D(listCoordinatesCartesianPointsOnInterface[0],listCoordinatesCartesianPointsOnInterface[1],listCoordinatesCartesianPointsOnInterface[2],listMagnitudeGradient)
    
 
    
def plotScatter2D(listPointsX,
                  listPointsY,
                  labelTitle,
                  labelX,
                  labelY):
                  
    plt.scatter(listPointsX, listPointsY)
    
    plt.title(labelTitle)
    
    plt.ylabel(labelX)
    
    plt.ylabel(labelY)
    
    plt.show()
    
    return        
    
#def plotMeltpoolSurface4D(pointsXref,pointsYref,pointsZref,valuesColorRef):
#    array_xy=[]
#    for ii in range(0,len(pointsXref)):
#        array_xy.append([])
#    for jj in range(0,len(pointsXref)):
#        array_xy[jj].append(pointsXref[jj])
#        array_xy[jj].append(pointsYref[jj])          
#
#    x_max=max(pointsXref)
#    x_min=min(pointsXref)
#    y_max=max(pointsYref)
#    y_min=min(pointsYref)
#    x = np.linspace(x_min,x_max,2000)
#    y = np.linspace(y_min,y_max,2000)    
#    grid_x,grid_y=np.meshgrid(x,y)
#    #listSamples = rand.sample(range(len(pointsYref)),min([len(pointsYref),2000])) 
#    #pointsXref_=[]
#    #pointsYref_=[]
#    #for ii in range(0,len(listSamples)):
#    #    pointsXref_.append(pointsXref[listSamples[ii]])
#    #    pointsYref_.append(pointsYref[listSamples[ii]])
#    #    
#    #grid_x,grid_y=np.meshgrid(pointsXref_,pointsYref_)
#    grid_z=griddata(np.array(array_xy),pointsZref,(grid_x,grid_y),method='nearest') 
#    grid_color=griddata(np.array(array_xy),valuesColorRef,(grid_x,grid_y),method='nearest')
#    
#    # eliminate points that are outside of the melt pool
#    #points=np.array([pointsXref,pointsYref]).transpose()
#    #hull = ConvexHull(points)
#    #for iPointX in range(0,len(x)):
#    #    for iPointY in range(0,len(y)):
#    #    np.append(points,[[x[iPointX],y[iPointY]]])
#    #    hullOuter = ConvexHull(points)
#    #    for simplex in hullOuter.simplices:
#            
#    
#    
#    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
#    C = grid_color
#    scamap = plt.cm.ScalarMappable(cmap='jet')
#    fcolors = scamap.to_rgba(C)
#    ax.plot_surface(grid_x, grid_y, grid_z, facecolors=fcolors, cmap='jet')
#    ax.set_xlabel("X")
#    ax.set_ylabel("Y")
#    ax.set_zlabel("Z")
#    fig.colorbar(scamap)
#    
#    # plot the boundary tracing the top of the melt pool
#    coordinateZmeltpoolTop=max(pointsZref)
#    points=np.array([pointsXref,pointsYref]).transpose()
#    hull = ConvexHull(points)
#    for simplex in hull.simplices:
#        ax.plot3D(points[simplex, 0], points[simplex, 1],np.array([coordinateZmeltpoolTop,coordinateZmeltpoolTop]),'black')
#    plt.show()
#    
#    return        
    
 
    
if __name__ == "__main__":

    main()