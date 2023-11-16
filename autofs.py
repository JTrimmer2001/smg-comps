from findclumpScript import findclump as fc
import os
import shutil

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

directory = 'testfiles' #name of the directory with the image files

for filename in os.listdir(directory):
    filename = filename.lower() #ensures everything is lower case

    cube = os.path.join(directory,filename) #produces a full file directory for reference purposes

    parts = filename.split('_') #separates the important bits of the file name
    
    if parts[1] == 'dirtycube':
        output = str('catalogues/'+parts[0]+'/'+parts[2])

        if os.path.isdir(output):
            continue
       
        #os.makedirs('./'+output) #temporary code to test shutil.copy
        fc(cube,output)

        for pdf in os.listdir(): #Trawls working directory looking for pdfs
            if pdf.endswith('.pdf') == True:
                shutil.copy(pdf,output) #copies all pdfs (ie the histogram files) into the targets catalogue folder


