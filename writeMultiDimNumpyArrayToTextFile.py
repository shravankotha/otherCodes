import numpy as np


data = np.array([[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[11, 12, 13], [14, 15, 16], [17, 18, 19]], [[21, 22, 23], [24, 25, 26], [27, 28, 29]]], np.int32)

with open('test1.dat', 'w') as outfile:
    # I'm writing a header here just for the sake of readability
    # Any line starting with "#" will be ignored by numpy.loadtxt
    outfile.write('# Array shape: {0}\n'.format(data.shape))
    
    # Iterating through a ndimensional array produces slices along
    # the last axis. This is equivalent to data[i,:,:] in this case
    for data_slice in data:
        print(data_slice)
        print('\n')
        # The formatting string indicates that I'm writing out
        # the values in left-justified columns 7 characters in width
        # with 2 decimal places.  
        np.savetxt(outfile, data_slice, fmt='%d')

        # Writing out a break to indicate different slices...
        outfile.write('# New slice\n')




#np.save('test.npy', arrayX, allow_pickle=False, fix_imports=True)
#Y=np.load('test.npy')
#print('Y:',Y)
#print('Y[0,:]:',Y[0][1][1])
