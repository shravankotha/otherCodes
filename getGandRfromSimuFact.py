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
import interpolate2D
import partitionCube
import hexElement 
import multiplyMatrices
import parseUnvFile_SimuFact 
import meltpoolCalculationsFromFE
import plotMeltpoolScatter4D

def main():

    nPoints, offsetTol, tolerance, methodInterpolation, timeStart  = 50, 0.005, 1E-6, 'linear', time.time()
    
    nArguments = len(sys.argv)
    
    if nArguments != 7:
    
        print('Error: Five command line arguments are expected -- (1) unv file name including path      \
                                                                  (2) liquidus temperature (same units as nodal temp values)   \
                                                                  (3) Laser direction (X or Y or Z)     \
                                                                  (4) Laser speed                       \
                                                                  (5) Melt pool normal direction')
                                                                  
        return
        
    nameFile, temperatureLiquidus = sys.argv[1], float(sys.argv[2])
    
    directionScanning, speedLaser = sys.argv[3], float(sys.argv[4])   
    
    directionMeltPoolNormal = sys.argv[5]
    
    iDimensionScanning = {"X":0, "Y":1, "Z":2}
    
    # -------------------------------------- Parse the unv file and obtain mesh information
    
    listNodes, listCoordinatesX, listCoordinatesY, listCoordinatesZ, listConnectivity, listTemperatures = parseUnvFile(nameFile)
    
    listElementIDs = [ii + 1 for ii in range(0,len(listConnectivity))]
                                                            
    listElementsInMeltPool = findElementsWithinTheMeltPool(listElementIDs, listConnectivity, listTemperatures, unitsTemperature, temperatureLiquidus)
    
    # -------------------------------------- Find the solid-liquid interface
    
    listElementsWithInterface = findElementsContainingInterfaceSolidLiquid(listElementIDs, listConnectivity, listTemperatures, temperatureLiquidus)
    
    listTemperaturesPointsOnInterface, listCoordinatesCartesianPointsOnInterface, \
            listVectorGradient, listMagnitudeGradient, listDendriteGrowthRate, listGrowthRate, listCoolingRate = ([] for ii in range(7))
    
    listCoordinatesCartesianPointsOnInterface = [[],[],[]], [[],[],[]]
    
    # -------------------------------------- Loop over all the interface elements and compute the interface points, thermal gradient, growth rate and cooling rate
    
    for iElement in range(0, len(listElementsWithInterface), 1):
    
        print("Processing element: ", iElement, " of ", str(len(listElementsWithInterface)), " -- elapsed time : ", str(time.time()-timeStart))
        
        listNodesConnectedElement = listConnectivity[listElementsWithInterface[iElement] - 1]
        
        listTemperaturesElement = [listTemperatures[iNode - 1] for iNode in listNodesConnectedElement]
        
        listCoordinatesNodalElement = [[],[],[]]
                       
        listCoordinatesNodalElement[0] = [listCoordinatesX[listNodesConnectedElement[iNode] - 1] for iNode in range(0, len(listNodesConnectedElement)]
                           
        listCoordinatesNodalElement[1] = [listCoordinatesY[listNodesConnectedElement[iNode] - 1] for iNode in range(0, len(listNodesConnectedElement)]
                           
        listCoordinatesNodalElement[2] = [listCoordinatesZ[listNodesConnectedElement[iNode] - 1] for iNode in range(0, len(listNodesConnectedElement)]
        
        listNaturalCoodsInterfaceCubeNodes, temperatureNodesNewWithInterface = findCubesEnclosingSolidLiquidInterfaceInAnElement(listNodesConnectedElement, listTemperaturesElement, temperatureLiquidus)
        
        # -------------------------------------- evaluate the centroid of each cube, temperature at the centroid and gradient cooling rate etc. at the centroid
                
        for iCube in range(0, len(listNaturalCoodsInterfaceCubeNodes)):
        
            listNaturalCoodsInterfaceCubeNodes_ = listNaturalCoodsInterfaceCubeNodes[iCube]
            
            listCoordinatesCartesianAtCentroid, listCoordinatesNaturalAtCentroid, temperatureAtCentroidgetCentroidValuesInterfaceCube = getCentroidValuesInterfaceCube(listTemperaturesElement,listCoordinatesNodalElement,listNaturalCoodsInterfaceCubeNodes)
                
            listCoordinatesCartesianPointsOnInterface[iDimension].append(listCartesianCoordinatesAtCentroid[iDimension])
            
            listTemperaturesPointsOnInterface.append(listTemperaturesAtCentroid)
            
            vectorGradient, magnitudeGradient = computeGradientScalarField(listCoordinatesNodalElement, listTemperaturesElement, listCoordinatesNaturalAtCentroid)
            
            for iDimension in range(0, 3):            
                
                listVectorGradient[iDimension].append(vectorGradient[iDimension])
                
            listMagnitudeGradient.append(math.sqrt(magnitudeGradient))
            
            listDendriteGrowthRate.append(speedLaser*abs(listVectorGradient[iDimensionScanning[directionScanning]][len(listMagnitudeGradient) - 1])/max([abs(listVectorGradient[0][len(listMagnitudeGradient) - 1]),abs(listVectorGradient[1][len(listMagnitudeGradient) - 1]),abs(listVectorGradient[2][len(listMagnitudeGradient) - 1])]))   # from: https://www.sciencedirect.com/science/article/abs/pii/S1005030216300615
            
            listGrowthRate.append(speedLaser*abs(listVectorGradient[iDimensionScanning[directionScanning]][len(listMagnitudeGradient) - 1])/math.sqrt(magnitudeGradient))
            
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
    
    lowestG, moderateG, highestG = listMagnitudeGradient[listIndices[0]], listMagnitudeGradient[listIndices[int(len(listIndices)/2)]], listMagnitudeGradient[listIndices[len(listIndices) - 1]]
        
    RforLowestG, RforModerateG, RforHighestG = listGrowthRate[listIndices[0]], listGrowthRate[listIndices[int(len(listIndices)/2)]], listGrowthRate[listIndices[len(listIndices) - 1]]
        
    VdforLowestG, VdforModerateG, VdforHighestG = listDendriteGrowthRate[listIndices[0]], listDendriteGrowthRate[listIndices[int(len(listIndices)/2)]], listDendriteGrowthRate[listIndices[len(listIndices) - 1]]
        
    # -------------------------------------- 
    listIndices = np.argsort(np.array(listGrowthRate))
    
    lowestR, moderateR, highestR = listGrowthRate[listIndices[0]], listGrowthRate[listIndices[int(len(listIndices)/2)]], listGrowthRate[listIndices[len(listIndices) - 1]]    
    
    GforLowestR, GforModerateR, GforHighestR = listMagnitudeGradient[listIndices[0]], listMagnitudeGradient[listIndices[int(len(listIndices)/2)]], listMagnitudeGradient[listIndices[len(listIndices) - 1]]    
    
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
    
    if directionScanning == "X" and directionMeltPoolNormal == "Y" or directionScanning == "Y" and directionMeltPoolNormal == "X":
    
        directionTransverse = "Z"
        
    elif directionScanning == "X" and directionMeltPoolNormal == "Z" or directionScanning == "Z" and directionMeltPoolNormal == "X":
    
        directionTransverse = "Y"

    elif directionScanning == "Y" and directionMeltPoolNormal == "Z" or directionScanning == "Z" and directionMeltPoolNormal == "Y":
    
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
    
if __name__ == "__main__":

    main()