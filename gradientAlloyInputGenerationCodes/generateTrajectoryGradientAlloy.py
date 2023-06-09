import sys
import os
import numpy as np
import math

def main():
   numLayers = 96
   thicknessLayer = 1.2
   velocityScan = 10.5
   arrayRadiiTrajectories = [0.5] + [1.0] + [1.5] + [2.0] + [2.5] + [3.0] + [3.5] + [3.5] 
   nTrajectories = len(arrayRadiiTrajectories)
   arrayAngleResolutionTrajectoryInDegrees = nTrajectories*[5]    # NOTE: should be an integer divisor of 360
   arrayPowers = 4*[8500000] + 4*[7900000] + 4*[7360000] + 4*[6840000] + 4*[6320000] + 76*[5780000]  # array size = numLayers
   arrayDwellTimes = 4*[4.5] + 91*[8.5]     # array size = numLayers - 1   
   
   #print(nTrajectories,len(arrayAngleResolutionTrajectoryInDegrees),len(arrayPowers),len(arrayDwellTimes))
   time = 0
   outfile = open('scanPath.inp', 'w')   
   outfile.write('** (1)time (2-4)X,Y,Z-coord of location of first event (5) first field value of first event (2kW laser)\n')
   outfile_step = open('abaqusSteps.inp', 'w')
   for iLayer in range(0,numLayers):    # goes from 0 to 95
        print('Processing layer : ' + str(iLayer))
        power = arrayPowers[iLayer]
        coord_2 = (iLayer+1)*thicknessLayer
        coord_prev_2 = coord_2
        nStepsPrinting = 0
        nStepsCooling = 0
        for iTrajectory in range(0,nTrajectories):            
            angleResolutionTrajectoryInDegrees = arrayAngleResolutionTrajectoryInDegrees[iTrajectory]
            nTimeStepsPerTrajectory = int(360/angleResolutionTrajectoryInDegrees) + 1            
            coord_prev_0 = arrayRadiiTrajectories[iTrajectory]
            coord_prev_1 = 0
            trajectoryTime = 0
            outfile.write(str(time) + ',' + str(coord_prev_0) + ',' + str(coord_prev_1) + ',' + str(coord_prev_2) + ',' + str(power) + '\n')
            for iTimeStepTrajectory in range(1,nTimeStepsPerTrajectory):
                coord_0 = arrayRadiiTrajectories[iTrajectory]*math.cos(iTimeStepTrajectory*angleResolutionTrajectoryInDegrees*math.pi/180)
                coord_1 = arrayRadiiTrajectories[iTrajectory]*math.sin(iTimeStepTrajectory*angleResolutionTrajectoryInDegrees*math.pi/180)
                #print('angle,coords: ',iTimeStepTrajectory*angleResolutionTrajectoryInDegrees,coord_0,coord_1)
                distance = math.sqrt((coord_0-coord_prev_0)**2 + (coord_1-coord_prev_1)**2 + (coord_2-coord_prev_2)**2)
                incrementTime = distance/velocityScan
                time = time + incrementTime
                trajectoryTime = trajectoryTime + incrementTime
                if iTrajectory == nTrajectories-1 and iTimeStepTrajectory == nTimeStepsPerTrajectory-1:
                    outfile.write(str(time) + ',' + str(coord_0) + ',' + str(coord_1) + ',' + str(coord_2) + ',' + str(0) + '\n')
                else:
                    outfile.write(str(time) + ',' + str(coord_0) + ',' + str(coord_1) + ',' + str(coord_2) + ',' + str(power) + '\n')
                coord_prev_0 = coord_0
                coord_prev_1 = coord_1
                coord_prev_2 = coord_2                
            if iTrajectory < nTrajectories-1:
                interTrajectoryTime = (arrayRadiiTrajectories[iTrajectory+1]-arrayRadiiTrajectories[iTrajectory])/velocityScan
                time = time + interTrajectoryTime
                trajectoryTime = trajectoryTime + interTrajectoryTime
            nStepsPrinting = nStepsPrinting + 1
            writeStep_Printing(outfile_step,trajectoryTime/10,trajectoryTime,1E-10,trajectoryTime/10,iLayer+1,nStepsPrinting)            
        if iLayer < numLayers-1:
            dwellTime = arrayDwellTimes[iLayer]
            time = time + dwellTime
            nStepsCooling = nStepsCooling + 1
            writeStep_Cooling(outfile_step,dwellTime/10,dwellTime,1E-10,dwellTime/10,iLayer+1,nStepsCooling)
        
   outfile.close()
   
def writeStep_Printing(fileHandle,initialTime,finalTime,minTime,maxTime,layerNumber,stepNumber):
    fileHandle.write('**-------------------------------------------------------------------------------------' + '\n')
    fileHandle.write('*STEP,INC=80000, NAME=Printing_layer_' + str(layerNumber) + '_traj_' + str(stepNumber) + ', UNSYMM=NO, EXTRAPOLATION=NO' + '\n')
    fileHandle.write('*HEAT TRANSFER\n')
    fileHandle.write(' ' + str(initialTime) + ',' + str(finalTime) + ',' + str(minTime) + ',' + str(maxTime) + '\n')
    fileHandle.write('*RESTART, WRITE, OVERLAY' + '\n')
    fileHandle.write('*ACTIVATE ELEMENTS,ACTIVATION="LDED.Printing"' + '\n')
    fileHandle.write(' "ABQ_AM.MaterialDeposition"' + '\n')
    fileHandle.write('*DFLUX' + '\n')
    fileHandle.write(' Printsubstrate, MBFNU, , "ABQ_AM.EnergyInput"' + '\n')
    fileHandle.write('*Solution Technique, Type=Quasi-Newton, Reform Kernel=8' + '\n')
    fileHandle.write('*OUTPUT,FIELD,FREQ=1' + '\n')
    fileHandle.write('*NODE OUTPUT' + '\n')
    fileHandle.write(' NT' + '\n') 
    fileHandle.write('*END STEP' + '\n')
    
def writeStep_Cooling(fileHandle,initialTime,finalTime,minTime,maxTime,layerNumber,stepNumber):
    fileHandle.write('**------------------------------------------------' + '\n')
    fileHandle.write('*STEP, INC=80000, NAME=Cooling_layer_' + str(layerNumber) + '_step_' + str(stepNumber) + ', UNSYMM=NO, EXTRAPOLATION=NO' + '\n')
    fileHandle.write('*HEAT TRANSFER' + '\n')
    fileHandle.write(' ' + str(initialTime) + ',' + str(finalTime) + ',' + str(minTime) + ',' + str(maxTime) + '\n')
    fileHandle.write('*RESTART, WRITE, FREQUENCY=10' + '\n')
    fileHandle.write('*Solution Technique, Type=Quasi-Newton, Reform Kernel=8' + '\n')
    fileHandle.write('*OUTPUT,FIELD,FREQ=1' + '\n')
    fileHandle.write('*NODE OUTPUT' + '\n')
    fileHandle.write(' NT' + '\n')
    fileHandle.write('*END STEP' + '\n')   
   
   
if __name__ == "__main__":

    main()