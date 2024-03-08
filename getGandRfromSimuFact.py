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

def main():

    methodInterpolation = 'linear'
    
    nPoints = 50
    
    shouldSaveFigures = True
    
    figureFormat = '.png'
    
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
    
    # -------------------------------------- Find the solid-liquid interface
    
    listElementsWithInterface = findInterfaceSolidLiquid(listConnectivity,
                                                         listTemperatures,
                                                         unitsTemperature,
                                                         temperatureLiquidus)
                                                         
    listElementsInMeltPool = findElementsInTheMeltPool(listConnectivity,
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
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    
    xlim_min = min(coordsGridScanningDirection[0]) - 6*offsetTol*(max(coordsGridScanningDirection[0])-min(coordsGridScanningDirection[0]))
    xlim_max = max(coordsGridScanningDirection[0]) + 6*offsetTol*(max(coordsGridScanningDirection[0])-min(coordsGridScanningDirection[0]))

    ylim_min = min(coordsGridNormalDirection[0]) - 6*offsetTol*(max(coordsGridNormalDirection[0])-min(coordsGridNormalDirection[0]))
    ylim_max = max(coordsGridNormalDirection[0]) + 6*offsetTol*(max(coordsGridNormalDirection[0])-min(coordsGridNormalDirection[0]))

    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    figNameToSave = './meltPool_2D' + str(figureFormat)
    if shouldSaveFigures == False:
        plt.title("Melt pool center line")
        plt.show()
    else:
        plt.savefig(figNameToSave,bbox_inches = 'tight',dpi = 300)
 
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
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    figNameToSave = './temperature_2D' + str(figureFormat)
    if shouldSaveFigures == False:
        plt.title("Temperature (K) ")    
        plt.show()
    else:
        plt.savefig(figNameToSave,bbox_inches = 'tight',dpi = 300)
    
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
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    figNameToSave = './thermalGradient_2D' + str(figureFormat)
    if shouldSaveFigures == False:
        plt.title("Thermal gradient")    
        plt.show()
    else:
        plt.savefig(figNameToSave,bbox_inches = 'tight',dpi = 300)
    
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
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    figNameToSave = './growthRate_2D' + str(figureFormat)
    if shouldSaveFigures == False:
        plt.title("Growth rate")    
        plt.show()
    else:
        plt.savefig(figNameToSave,bbox_inches = 'tight',dpi = 300)
    
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
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    figNameToSave = './dendriteGrowthRate_2D' + str(figureFormat)
    if shouldSaveFigures == False:
        plt.title("Dendrite growth rate")    
        plt.show()
    else:
        plt.savefig(figNameToSave,bbox_inches = 'tight',dpi = 300)
    
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
    plt.xlabel("Length (m)")
    plt.ylabel("Depth (m)")
    plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
    figNameToSave = './coolingRate_2D' + str(figureFormat)
    if shouldSaveFigures == False:
        plt.title("Cooling rate")    
        plt.show()
    else:
        plt.savefig(figNameToSave,bbox_inches = 'tight',dpi = 300)
    
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
                          './meltPool_3D' + str(figureFormat),
                          shouldSaveFigures = shouldSaveFigures,                          
                          shouldTurnOffGrid = False,
                          shouldTurnOffAxes = False,
                        )
                          
    plotMeltpoolScatter4D(listCoordinatesCartesianPointsOnInterface[0],
                          listCoordinatesCartesianPointsOnInterface[1],
                          listCoordinatesCartesianPointsOnInterface[2],
                          listMagnitudeGradient,
                          directionMeltPoolNormal,
                          "Thermal Gradient (G)",
                          "X (m)",
                          "Y (m)",
                          "Z (m)",
                          './thermalGradient_3D' + str(figureFormat),
                          shouldSaveFigures = shouldSaveFigures,
                          shouldTurnOffGrid = True,
                          shouldTurnOffAxes = True
                          )
                          
    plotMeltpoolScatter4D(listCoordinatesCartesianPointsOnInterface[0],
                          listCoordinatesCartesianPointsOnInterface[1],
                          listCoordinatesCartesianPointsOnInterface[2],
                          listGrowthRate,
                          directionMeltPoolNormal,
                          "Growth Rate (R)",
                          "X (m)",
                          "Y (m)",
                          "Z (m)",
                          './growthRate_3D' + str(figureFormat),
                          shouldSaveFigures = shouldSaveFigures,
                          shouldTurnOffGrid = True,
                          shouldTurnOffAxes = True
                          )

    plotMeltpoolScatter4D(listCoordinatesCartesianPointsOnInterface[0],
                          listCoordinatesCartesianPointsOnInterface[1],
                          listCoordinatesCartesianPointsOnInterface[2],
                          listDendriteGrowthRate,
                          directionMeltPoolNormal,
                          "Dendrite Growth Rate (R)",
                          "X (m)",
                          "Y (m)",
                          "Z (m)",
                          './dendriteGrowthRate_3D' + str(figureFormat),
                          shouldSaveFigures = shouldSaveFigures,
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
                          './coolingRate_3D' + str(figureFormat),
                          shouldSaveFigures = shouldSaveFigures,
                          shouldTurnOffGrid = True,
                          shouldTurnOffAxes = True
                          )

    # scatter plot the G vs R
    #plotScatter2D(listMagnitudeGradient,listGrowthRate,"G versus R","G (K/m)","R (m/s)")
    #plotMeltpoolSurface4D(listCoordinatesCartesianPointsOnInterface[0],listCoordinatesCartesianPointsOnInterface[1],listCoordinatesCartesianPointsOnInterface[2],listMagnitudeGradient)
    
    
def interpolate2D(listPointsXref, listPointsYref, listPointsZref, \
                xMinToBeInterpolated, xMaxToBeInterpolated, \
                yMinToBeInterpolated, yMaxToBeinterpolated, \
                nPoints, methodInterpolation):
                
    array_xy = []
    
    for ii in range(0, len(listPointsXref)):
    
        array_xy.append([])
        
    for jj in range(0, len(listPointsXref)):
    
        array_xy[jj].append(listPointsXref[jj])
        
        array_xy[jj].append(listPointsYref[jj])                 

    x = np.linspace(xMinToBeInterpolated, xMaxToBeInterpolated, nPoints)
    
    y = np.linspace(yMinToBeInterpolated, yMaxToBeinterpolated, nPoints)    
    
    grid_x, grid_y = np.meshgrid(x, y)
    
    grid_z = griddata(np.array(array_xy), listPointsZref, (grid_x, grid_y), method = methodInterpolation) 
    
    return grid_x, grid_y, grid_z 
    
    
def plotScatter2D(listPointsX,
                  listPointsY,
                  labelTitle,
                  labelX,
                  labelY,
                  pathFigNameToSave,
                  shouldSaveFigures):
                  
    plt.scatter(listPointsX, listPointsY)    
    
    plt.ylabel(labelX)
    
    plt.ylabel(labelY)
    
    if shouldSaveFigures == False:
        plt.title(labelTitle)
        plt.show()
    else:
        plt.savefig(pathFigNameToSave,bbox_inches = 'tight',dpi = 300)
    
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
    
    
    
def plotMeltpoolScatter4D(listPointsX,
                          listPointsY,
                          listPointsZ,
                          listPointsColor,
                          directionMeltPoolNormal,
                          labelTitle,
                          labelX,
                          labelY,
                          labelZ,
                          pathFigNameToSave,
                          shouldSaveFigures,
                          shouldTurnOffGrid,
                          shouldTurnOffAxes
                          ):
                          
    plt.rcParams.update({'font.family':'Times New Roman',
                         'font.weight': 'light',
                         'font.size': 14})
                         
    plt.rcParams.update()
    
    fig = plt.figure()
    
    ax = plt.axes(projection = '3d')
    
    img = ax.scatter(listPointsX, listPointsY, listPointsZ, c = listPointsColor, cmap = 'jet', zorder = -1)
    
    fig.colorbar(img)        
    
    ax.set_xlabel(labelX, labelpad = 20, fontweight = 'normal')
    
    ax.set_ylabel(labelY, labelpad = 20, fontweight = 'normal')
    
    ax.set_zlabel(labelZ, labelpad = 20, fontweight = 'normal')
    
    # plot the boundary tracing the top of the melt pool
    
    if directionMeltPoolNormal == "X":
    
        coordinateMeltpoolTop = max(listPointsX)
        
        points = np.array([listPointsY, listPointsZ]).transpose()
        
    elif directionMeltPoolNormal == "Y":
    
        coordinateMeltpoolTop = max(listPointsY)
        
        points = np.array([listPointsX, listPointsZ]).transpose()

    elif directionMeltPoolNormal == "Z":
    
        coordinateMeltpoolTop = max(listPointsZ)
        
        points = np.array([listPointsX, listPointsY]).transpose()
    
    hull = ConvexHull(points)
    
    for simplex in hull.simplices:
    
        ax.plot3D(points[simplex, 0], points[simplex, 1], np.array([coordinateMeltpoolTop, coordinateMeltpoolTop]), 'black', zorder = 5)
    
    ax.view_init(elev = 50, azim = 205)
    
    if shouldTurnOffGrid == True:
    
        plt.grid(False)
        
    if shouldTurnOffAxes == True:
    
        plt.axis('off')
        
    if shouldSaveFigures == False: 
        ax.set_title(labelTitle, fontsize = 18, fontweight = 'bold')    
        plt.show()
    else:
        plt.savefig(pathFigNameToSave,bbox_inches = 'tight',dpi = 300)
    
    return

def findInterfacesSolidLiquidInAnElement(listNodesConnectedElement,
                                         listTemperaturesElement,
                                         listCoordinatesNodalElement,
                                         temperatureLiquidus):
                                         
    listXi, listEta, listZeta = coordsNaturalElementHex()
    
    coordinatesNodesOld = []
    
    coordinatesNodesOld.append([listXi, listEta, listZeta])
    
    nDivisions = 8
    
    maxPartitionedCubesToKeep = 1
    
    for iDivision in range(0, nDivisions):
    
        coordinatesNodesNewWithInterface = []
        
        temperatureNodesNewWithInterface = []
        
        for iCube in range(0, len(coordinatesNodesOld)):
        
            coordinatesNodesNew, coordinatesNewElementsHex = partitionCube(coordinatesNodesOld[iCube])
            
            for iCubeNew in range(0, len(coordinatesNewElementsHex)):
            
                coordinatesNodeElemental = coordinatesNewElementsHex[iCubeNew]
                
                isFoundNodeWithLowerThanLiquidus = False
                
                isFoundNodeWithHigherThanLiquidus = False
                
                listTemperatures = []
                
                for iNode in range(0, len(coordinatesNodeElemental[0])):
                
                    coordinatesNatural = [coordinatesNodeElemental[0][iNode], coordinatesNodeElemental[1][iNode], coordinatesNodeElemental[2][iNode]]
                    
                    listTemperatures.append(getQuantityInsideTheElement(listTemperaturesElement, coordinatesNatural))
                    
                    if listTemperatures[iNode] < temperatureLiquidus:
                    
                        isFoundNodeWithLowerThanLiquidus = True
                        
                    elif listTemperatures[iNode] >= temperatureLiquidus:
                    
                        isFoundNodeWithHigherThanLiquidus = True
                    
                if isFoundNodeWithLowerThanLiquidus == True and isFoundNodeWithHigherThanLiquidus == True:
                
                    coordinatesNodesNewWithInterface.append(coordinatesNodeElemental)
                    
                    temperatureNodesNewWithInterface.append(listTemperatures)
                    
        listCubeIDs = rand.sample(range(len(coordinatesNodesNewWithInterface)), min([len(coordinatesNodesNewWithInterface), maxPartitionedCubesToKeep])) 
        
        coordinatesNodesOld = []
        
        for ii in range(0, len(listCubeIDs)):
        
            coordinatesNodesOld.append(coordinatesNodesNewWithInterface[listCubeIDs[ii]])

    return coordinatesNodesNewWithInterface, temperatureNodesNewWithInterface
    
def getQuantityInsideTheElement(quantitiesNodal, coordinatesNatural):

    shapeFunctions = evaluateShapeFunctions(coordinatesNatural)
    
    quantity = 0
    
    for iNode in range(0,len(quantitiesNodal)):
    
        quantity = quantity + shapeFunctions[iNode]*quantitiesNodal[iNode]
        
    return quantity
    
def evaluateShapeFunctions(coordinatesNatural):

    listXi, listEta, listZeta = coordsNaturalElementHex()
    
    xi = coordinatesNatural[0]
    
    eta = coordinatesNatural[1]
    
    zeta = coordinatesNatural[2]
    
    shapeFunctions = []
    
    for iNode in range(0, 8):
    
        shapeFunctions.append((1/8)*(1 + listXi[iNode]*xi)*(1 + listEta[iNode]*eta)*(1 + listZeta[iNode]*zeta))
        
    return shapeFunctions
        
    
def partitionCube(coordinatesNodesOld):

    coordinatesNodesNew = []
    
    for iDimension in range(0, 3):
    
        coordinatesNodesNew.append([])
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][0] + coordinatesNodesOld[iDimension][1]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][1] + coordinatesNodesOld[iDimension][2]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][2] + coordinatesNodesOld[iDimension][3]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][3] + coordinatesNodesOld[iDimension][0]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][4] + coordinatesNodesOld[iDimension][5]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][5] + coordinatesNodesOld[iDimension][6]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][6] + coordinatesNodesOld[iDimension][7]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][7] + coordinatesNodesOld[iDimension][4]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][1] + coordinatesNodesOld[iDimension][5]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][2] + coordinatesNodesOld[iDimension][6]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][3] + coordinatesNodesOld[iDimension][7]))
        
        coordinatesNodesNew[iDimension].append((1/2)*(coordinatesNodesOld[iDimension][0] + coordinatesNodesOld[iDimension][4]))
        
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][0] + coordinatesNodesOld[iDimension][1] + coordinatesNodesOld[iDimension][2] + coordinatesNodesOld[iDimension][3]))
        
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][4] + coordinatesNodesOld[iDimension][5] + coordinatesNodesOld[iDimension][6] + coordinatesNodesOld[iDimension][7]))
        
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][0] + coordinatesNodesOld[iDimension][4] + coordinatesNodesOld[iDimension][5] + coordinatesNodesOld[iDimension][1]))
        
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][3] + coordinatesNodesOld[iDimension][7] + coordinatesNodesOld[iDimension][6] + coordinatesNodesOld[iDimension][2]))
        
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][0] + coordinatesNodesOld[iDimension][3] + coordinatesNodesOld[iDimension][7] + coordinatesNodesOld[iDimension][4]))
        
        coordinatesNodesNew[iDimension].append((1/4)*(coordinatesNodesOld[iDimension][1] + coordinatesNodesOld[iDimension][2] + coordinatesNodesOld[iDimension][6] + coordinatesNodesOld[iDimension][5]))
        
        coordinatesNodesNew[iDimension].append((1/8)*(coordinatesNodesOld[iDimension][0] + coordinatesNodesOld[iDimension][1] + coordinatesNodesOld[iDimension][2] + coordinatesNodesOld[iDimension][3] +
                                                  coordinatesNodesOld[iDimension][4] + coordinatesNodesOld[iDimension][5] + coordinatesNodesOld[iDimension][6] + coordinatesNodesOld[iDimension][7]))
    
    coordinatesNewElementsHex = []
    
    for iCube in range(0, 8):
    
        coordinatesNewElementsHex.append([])
        
        for iDimension in range(0, 3):
        
            coordinatesNewElementsHex[iCube].append([])
        
    for iDimension in range(0, 3):
    
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
    
def dShapeFunction_dCaterianCoords(coordinatesCartesianNodal, coordinatesNaturalPoint):

    Jacobian, dN_dXi, dN_dEta, dN_dZeta = getJacobianHexElement(coordinatesCartesianNodal, coordinatesNaturalPoint)
    
    matJacobian = np.matrix(Jacobian)
    
    matInvJacobian = matJacobian.I
    
    invJacobian = matInvJacobian.tolist()
    
    dN_dxyz = []
    
    for iNode in range(0, 8):
    
        dN_dxyz.append(np.matmul(matInvJacobian,np.matrix([dN_dXi[iNode], dN_dEta[iNode], dN_dZeta[iNode]]).transpose()).transpose().tolist()[0])
        
    return dN_dxyz    
       
def getJacobianHexElement(coordinatesCartesianNodal, coordinatesNaturalPoint):

    x, y, z = ([] for ii in range(3))
    
    for iNode in range(0,8):
    
        x.append(coordinatesCartesianNodal[0][iNode])
        
        y.append(coordinatesCartesianNodal[1][iNode])
        
        z.append(coordinatesCartesianNodal[2][iNode])
        
    dN_dXi, dN_dEta, dN_dZeta = dShapeFunction_dNaturalCoords(coordinatesNaturalPoint)

    matrixA = []
    
    matrixB = []
    
    for ii in range(0, 8):
    
        matrixB.append([])
        
    for ii in range(0, 3):
    
        matrixA.append([])
        
        for jj in range(0, 8):
        
            if ii == 0:
            
                matrixA[ii].append(dN_dXi[jj])
                
                matrixB[jj].append(x[jj])
                
            elif ii == 1:
            
                matrixA[ii].append(dN_dEta[jj])
                
                matrixB[jj].append(y[jj])
                
            elif ii == 2:
            
                matrixA[ii].append(dN_dZeta[jj])
                
                matrixB[jj].append(z[jj])

    Jacobian = multiplyMatrices(matrixA, matrixB)
    
    return Jacobian, dN_dXi, dN_dEta, dN_dZeta
    
def dShapeFunction_dNaturalCoords(coordinatesNatural):

    xi = coordinatesNatural[0]
    
    eta = coordinatesNatural[1]
    
    zeta = coordinatesNatural[2]
    
    listXi, listEta, listZeta = coordsNaturalElementHex()
    
    dN_dXi, dN_dEta, dN_dZeta=([] for ii in range(3))
    
    for iNode in range(0, 8):
    
        dN_dXi.append((1/8)*listXi[iNode]*(1 + listEta[iNode]*eta)*(1 + listZeta[iNode]*zeta))
        
        dN_dEta.append((1/8)*(1 + listXi[iNode]*xi)*listEta[iNode]*(1 + listZeta[iNode]*zeta))
        
        dN_dZeta.append((1/8)*(1 + listXi[iNode]*xi)*(1 + listEta[iNode]*eta)*listZeta[iNode])
    
    return dN_dXi, dN_dEta, dN_dZeta   
  

def coordsNaturalElementHex():

    listXi, listEta, listZeta=([] for ii in range(3))
    
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



def findInterfaceSolidLiquid(listConnectivity,
                             listTemperatures,
                             unitsTemperature,
                             temperatureLiquidus):
                             
    listElementsWithInterface = []
    
    for iElement in range(0, len(listConnectivity)):
    
        listConnectedNodes = listConnectivity[iElement]
        
        isFoundNodeWithLowerThanLiquidus = False
        
        isFoundNodeWithHigherThanLiquidus = False
        
        for iNode in range(0, len(listConnectedNodes)):
        
            idNode = listConnectedNodes[iNode] - 1
            
            temperature = listTemperatures[idNode]
            
            if temperature < temperatureLiquidus:
            
                isFoundNodeWithLowerThanLiquidus = True
                
                idNodeWithLowerThanLiquidus = idNode
                
            elif temperature >= temperatureLiquidus:
            
                isFoundNodeWithHigherThanLiquidus = True
                
                idNodeWithHigherThanLiquidus = idNode
                
            else:
            
                raise ValueError("Error: Something wrong with the temperature value")
                
                
            if isFoundNodeWithLowerThanLiquidus == True and isFoundNodeWithHigherThanLiquidus == True:
            
                listElementsWithInterface.append(iElement + 1)
                
                break
                
    return listElementsWithInterface
    
def findElementsInTheMeltPool(listConnectivity,
                             listTemperatures,
                             unitsTemperature,
                             temperatureLiquidus):
                             
    listElementsInMeltPool = []
    
    for iElement in range(0, len(listConnectivity)):
    
        listConnectedNodes = listConnectivity[iElement]
        
        isElementInsideTheMeltPool = True
        
        for iNode in range(0, len(listConnectedNodes)):
        
            idNode = listConnectedNodes[iNode] - 1
            
            temperature = listTemperatures[idNode]
            
            if temperature < temperatureLiquidus:
            
                isElementInsideTheMeltPool = False
            
                break                
                                
        if isElementInsideTheMeltPool == True:
            
            listElementsInMeltPool.append(iElement + 1)                
                
    return listElementsInMeltPool
                
def multiplyMatrices(A, B):

    if len(A[0]) != len(B):
    
        raise ValuError("Error: The dimensions of the matrices are inconsistent for multiplication")
        
    C = []
    
    for ii in range(0, len(A)):
    
        C.append([])
        
        C[ii] = [0]*len(B[0])
        
        for jj in range(0, len(B[0])):       
        
            for kk in range(0, len(A[ii])):
            
                C[ii][jj] = C[ii][jj] + A[ii][kk]*B[kk][jj]
                
    return C
        
            
    
def parseUnvFile(nameFile):

    listNodesTemperature, listNodesNodalCoordinates, listIdElementsConnectivity = ([] for ii in range(3))
    
    with open(nameFile) as fileCurrent:
    
        dataTemp = fileCurrent.readlines()
        
        isFoundKeywordNodes = False
        
        isFoundKeywordConnectivity = False
        
        isFoundKeywordTemp = False
        
        shouldReadNodes = False
        
        shouldKeepReadingNodes = True
        
        shouldKeepReadingConnectivity = True
        
        count = 0
        
        for iData in range(0, len(dataTemp)):
        
            # Reading the nodal coordinates
            
            if str(dataTemp[iData].strip()) == "Contained bodies:":
            
                isFoundKeywordNodes = True
                
                count = 0
                
            if isFoundKeywordNodes == True and shouldKeepReadingNodes == True:
            
                count = count + 1
                
                if count >= 6:
                
                    listNodesNodalCoordinates.append(dataTemp[iData].strip().split())
                    
                if str(dataTemp[iData].strip()) == "-1" and str(dataTemp[iData+1].strip()) == "-1":
                
                    shouldKeepReadingNodes = False
                
            # Reading element connectivity
            
            if str(dataTemp[iData].strip()) == "2412" and str(dataTemp[iData-1].strip()) == "-1":
            
                isFoundKeywordConnectivity = True
                
                count = 0
                
            if isFoundKeywordConnectivity == True and shouldKeepReadingConnectivity == True:
            
                count = count + 1
                
                if count > 1:
                
                    listIdElementsConnectivity.append(dataTemp[iData].strip().split())
                    
                if str(dataTemp[iData].strip()) == "-1" and str(dataTemp[iData+1].strip()) == "-1":
                
                    shouldKeepReadingConnectivity = False
                    
            # Reading the temperature data
            
            if str(dataTemp[iData].strip()) == "Temperature":
            
                isFoundKeywordTemp = True
                
                count = 0
                
            if isFoundKeywordTemp == True:
            
                if count <= 9:    # skipping 11 lines after Temperature keyword is found
                
                    count = count + 1
                    
                    continue
                    
                else:
                
                    if iData != len(dataTemp) - 1:
                    
                        listNodesTemperature.append(dataTemp[iData].strip())
                        
        # Do some checks and write coordinates and temperatures into lists
        
        if isFoundKeywordNodes == False or isFoundKeywordConnectivity == False or isFoundKeywordTemp == False:
        
            sys.exit("Error: Nodal coords or connectivity or Temperature keyword/information is not found in the file")
            
        else:
        
            listNodes, listCoordinatesX, listCoordinatesY, listCoordinatesZ, listElements, listConnectivity = ([] for ii in range(6))
            
            # Place nodal coordinates into arrays
            
            for i in range(0, len(listNodesNodalCoordinates) - 1, 2):
            
                listNodes.append(int(listNodesNodalCoordinates[i][0]))
                
            for i in range(1, len(listNodesNodalCoordinates) - 1, 2):
            
                listCoordinatesX.append(float(listNodesNodalCoordinates[i][0]))
                
                listCoordinatesY.append(float(listNodesNodalCoordinates[i][1]))
                
                listCoordinatesZ.append(float(listNodesNodalCoordinates[i][2]))
                
            # Place element connectivity into arrays
            
            for i in range(0, len(listIdElementsConnectivity) - 1, 2):
            
                listElements.append(int(listIdElementsConnectivity[i][0]))
                
            count = 0
            
            for i in range(1, len(listIdElementsConnectivity) - 1, 2):
            
                count = count + 1
                
                listConnectivity.append([])
                
                for j in range(0, 8):
                
                    listConnectivity[count - 1].append(int(listIdElementsConnectivity[i][j]))
                    
            # Place temperature into arrays
            
            listNodesTempData = [int(i) for i in listNodesTemperature[0:len(listNodesTemperature):2]]
            
            listTemperatures = [float(i) for i in listNodesTemperature[1:len(listNodesTemperature):2]]
            
    return listNodes, listCoordinatesX, listCoordinatesY, listCoordinatesZ, listConnectivity, listTemperatures
    
if __name__ == "__main__":

    main()