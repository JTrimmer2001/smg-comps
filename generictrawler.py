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

    files = []

    if form != 'none':
        if form.startswith('.') == False:
            print('Value error: . missing from beginning of format so this will be added in post. Please add a . in future runs.')
        
            format = str('.'+form)

    for filename in os.listdir(path):

        if form == 'none' and kword == 'none':

            files.append(filename)
        
        elif form != 'none' and filename.endswith(form) == True:
            if kword != 'none' and filename.count(kword)>0:
                files.append(filename)
            
            elif kword == 'none':
                files.append(filename)

        elif kword != 'none' and filename.count(kword) >0:
            if form != 'none' and filename.endswith(form) == True:
                files.append(filename)
            
            elif form == 'none':
                files.append(filename)

    if len(files) == 0:
        print('No entries found on this run, perhaps your search terms are wrong/too strict?')
        files = str('none')

    elif len(files) == 1:
        files=str(files[0])
        
    return(files)

def filefilter():
    dir = trawler('catalogues')

    for i in dir:
        cats = trawler('catalogues/'+i,kword='_wfidelity')

        for m in cats:
            su.copy('catalogues/'+i+'/'+m,'mastercats/'+i+m)
