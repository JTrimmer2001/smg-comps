import shutil as su
import os


def trawler(path, form = 'none', kword = 'none'):
    '''
    Script to trawl a directory for files, returns a list of the file names

    Inputs-
        path (str): path to directory 
        format (str): (default: none) specific type of file to look for, ie .pdf or .doc
        kword (str): (default: none) only files containing the keyword will be recorded

    Outputs-
        files: list of qualifying files found within the directory
    '''

    filename=str(filename)
    form = str(form)
    kword = str(kword)
    files = []

    if form != 'none':
        if form.startswith('.') == False:
            print('Value error: . missing from beginning of format so this will be added in post. Please add a . in future runs.')
        
            format = str('.'+form)

    for filename in os.listdir(path):

        if form == 'none' and kword == 'none':

            files.append(filename)
        
        elif form != 'none' and filename.endswith(form) == True:
            if kword != 'none' and filename.count(kword)==0:
                files.append(filename)
            
            elif kword == 'none':
                files.append(filename)

        elif kword != 'none' and filename.count(kword) == 0:
            if form != 'none' and filename.endswith(form) == True:
                files.append(filename)
            
            elif form == 'none':
                files.append(filename)