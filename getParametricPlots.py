import sys
import os
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
import numpy as np
import math



idColumnX = 1
idColumnY = 3
idColumnData = 8        # 6-thermalGradient, 7-growthRate, 8-coolingRate
multiplierForXaxisPlot = 1E6
multiplierForYaxisPlot = 1E6
sizeFont1 = 28
sizeFont2 = 24
scatterPlotMarkerSize = 50
fontTimes = {'fontname':'Times New Roman'}
sizeFont = sizeFont1
plt.rcParams.update({'font.size': sizeFont2})


if idColumnData == 6:
    labelPlot = 'Thermal Gradient (K/m)'
elif idColumnData == 7:
    labelPlot = 'Growth Rate (m/s)'
elif idColumnData == 8:
    labelPlot = 'Cooling Rate (K/s)'
else:
    raise ValueError('User Error: This column id is not valid')

coordsGridScanningDirection_All = []
coordsGridNormalDirection_All = []
gridData_All = []

#
out_paths = ['G:\My Drive\AFOSR project\C103_simulations_presentations\Data for C103 G.R plots - old\\125W-250mm-s-Data for G.R\\' + 'coords_depth_Temp_G_R_coolingRate_alongCenterline.out', \
'G:\My Drive\AFOSR project\C103_simulations_presentations\Data for C103 G.R plots - old\\175W-600mm-s-Data for G.R\\' + 'coords_depth_Temp_G_R_coolingRate_alongCenterline.out', \
'G:\My Drive\AFOSR project\C103_simulations_presentations\Data for C103 G.R plots - old\\230W-1200mm-s-Data for G.R\\' + 'coords_depth_Temp_G_R_coolingRate_alongCenterline.out']


for iFile in range(0,len(out_paths)):
    out_path = out_paths[iFile]
    file = open(out_path,mode='r') 
    lines = file.readlines()[1:]    # skip the header

    for line in lines:  
        listStrings = line.strip().split(' ')
        listData = []
        for iString in range(0,len(listStrings)):
            if listStrings[iString] != '':
                listData.append(float(listStrings[iString]))
        coordsGridScanningDirection_All.append(multiplierForXaxisPlot*float(listData[idColumnX-1]))
        coordsGridNormalDirection_All.append(multiplierForYaxisPlot*(listData[idColumnY-1]))
        gridData_All.append(listData[idColumnData-1])
    file.close()    

print(gridData_All)

plt.figure()

#plt.rc('xtick', labelsize=25)
#plt.rc('ytick', labelsize=25)
sc = plt.scatter(coordsGridScanningDirection_All, coordsGridNormalDirection_All, c = gridData_All, cmap = "jet", s=scatterPlotMarkerSize)
plt.colorbar()
#plt.title(labelPlot,**fontTimes)
plt.xlabel("Length ($\mu$m)",**fontTimes,fontsize=sizeFont)
plt.ylabel("Depth ($\mu$m)",**fontTimes,fontsize=sizeFont)
plt.tight_layout()

#ax.set_aspect('equal', adjustable='box')

#plt.axis([xlim_min, xlim_max, ylim_min, ylim_max])
plt.show()