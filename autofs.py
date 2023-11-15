from findclumpScript import findclump as fc
import os

# directory = 'testfiles'

# for filename in os.listdir(directory):
#     f = os.path.join(directory,filename)

#     if os.path.isfile(f):
#         print(f)

# Above is the example for iterating through filenames graciously provided by geeks4geeks
########################################################################################

'''
So the plan is:
    > Find the next file in the directory
    > Recover the source ID (alessXX)
    > Recover whether the file is a cube or not
    > Recover the spectral window (spwXXXX or allspw)
    > Form a folder address from this (concatenate??)
    > Run the findclumps script for this output and input

    Cant be too hard right?
'''

directory = 'testfiles'

for filename in os.listdir(directory):
    filename = filename.lower()

    cube = os.path.join(directory,filename)

    parts = filename.split('_')
    
    if parts[1] == 'dirtycube':
        output = str('catalogues/'+parts[0]+'/'+parts[2])
        print(output)
        fc(cube,output)
    else:
        print('flat :(')
