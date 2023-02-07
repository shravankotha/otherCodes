import sys
import os
import numpy as np
import math

def main():
    nArguments = len(sys.argv)
    if nArguments != 4:
    
        print('Error: Three command line arguments are expected -- (1) Material name\
        (2) Angle in degrees from SD direction (in SD-BD plane)\
        (3) Temperature (in C)')
        return
        
    nameMaterial = sys.argv[1]
    angleDegrees = float(sys.argv[2])
    temperature = float(sys.argv[3])
    
    if nameMaterial == '316L':
        if temperature < 23 or temperature > 400:
            print('Input Error: temperature must be between 23C and 400C for this material')
            return             
        coefficientF = 0.002*temperature + 0.7135
        coefficientG = 0.002*temperature + 0.7135
        coefficientH = -0.002*temperature + 0.2865
        coefficientL = 1.60
        coefficientM = 1.58
        coefficientN = 1.60
        stressYieldSD = -0.366*temperature + 481.4        
    elif nameMaterial == '17-4PH':
        if temperature < 23 or temperature > 1000:
            print('Input Error: temperature must be between 23C and 1000C for this material')
            return     
        coefficientF = -(4E-5)*temperature + 0.4462
        coefficientG = -(4E-5)*temperature + 0.4462
        coefficientH = 0.014*temperature + 0.535
        coefficientL = 1.60
        coefficientM = 1.58
        coefficientN = 1.60
        stressYieldSD = -0.8831*temperature + 990.15
    else:
        print('Input Error: The given material database is not availabe. Available materials - 316L and 17-4PH')
        return
    
    mappingDirectionsAnisotropy = {"SD":1, "GD":2, "BD":3}
    
    # -------------------------------------- Transform the unknown uni-axial stress along theta to anisotropy frame
    angleRadians = math.radians(angleDegrees)
    R11 = (np.cos(angleRadians))**2
    R12 = (np.sin(angleRadians))**2
    R13 = np.sin(2*angleRadians)
    R21 = R12
    R22 = R11
    R23 = -R13
    R31 = -np.sin(angleRadians)*np.cos(angleRadians)
    R32 = -R31
    R33 = np.cos(2*angleRadians)
    
    matrixR = np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]])
    matrixInvR = np.linalg.inv(matrixR)
    
    stressUniaxialYield = stressYieldSD/math.sqrt(coefficientF*(matrixInvR[1][0])**2 + coefficientG*(matrixInvR[1][0] - \
                            matrixInvR[0][0])**2 + coefficientH*(matrixInvR[0][0])**2 + 2*coefficientM*(matrixInvR[2][0])**2)
    print("Computed yield stress: ", stressUniaxialYield, " MPa")
    
if __name__ == "__main__":

    main()