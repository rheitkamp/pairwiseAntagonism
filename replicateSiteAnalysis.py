import os

this_dir = os.getcwd()
os.chdir(this_dir)

##Read pairwise interactions
def loadPairwiseAntagonism(inputFilename): 
    '''
    Read in a dictionary saved as a pairwise interaction file.

    Lines:
        (Plated Isolate)    (Pin Isolate)   (Interaction score)
    '''
    with open(inputFilename, 'r') as inputFile:
        lines = inputFile.readlines()
    
    dictionary = {}
    
    for line in lines:
        line = line.strip()
        parseLine = line.split('\t')
        
        if parseLine[2] == 'True':
            parseLine[2] = True
        elif parseLine[2] == 'False':
            parseLine[2] = False

        try:
            dictionary[parseLine[0]][parseLine[1]] = parseLine[2]
        except:
            dictionary[parseLine[0]] = {}
            dictionary[parseLine[0]][parseLine[1]] = parseLine[2]
    return dictionary

def dictionaryToMatrix(outputFilename, dictionary):
    '''
    '''    
    report = dictionary     

    columns = []
    isolates = report.keys()
    for isolate2 in report[isolates[0]]:
            if isolate2[-1] != '*':
                columns.append(isolate2)
    columns.sort()
    #print columns

    lines = []
    
    firstRow = ''
    for isolate2 in columns:
        firstRow += '\t' + str(isolate2)
    firstRow += '\n'
    lines.append(firstRow)
    
    isolates.sort()
    
    for isolate in isolates:
        row = str(isolate)
            
        for isolate2 in columns:
            #print isolate, killer
            score = report[isolate][isolate2]
            #print interaction
            row += '\t' + str(int(score)) 
        row += '\n'
        lines.append(row)
        
    with open(outputFilename, 'w') as outputFile:
        outputFile.writelines(lines)


filesDir = 'E:\\Antagonism\\First 94 Acinetobacter Antagonism Results\\Antagonism Data Sheets\\Analysis\\SampleData\\AcinetobacterData\\replicate-sites'
os.chdir(filesDir)
files = os.listdir(os.getcwd())
outfiles = ('rep1L.csv','rep2L.csv','rep3L.csv','rep1R.csv','rep2R.csv','rep3R.csv')
index = 0
for eachFile in files:
    siteRep = loadPairwiseAntagonism(eachFile)
    os.chdir(this_dir)
    dictionaryToMatrix(outfiles[index], siteRep)
    index+=1
    os.chdir(filesDir)
