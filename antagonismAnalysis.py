import os

#os.chdir('E:\\Antagonism\\First 94 Acinetobacter Antagonism Results\\Antagonism Data Sheets\\Analysis\\posterOutput') #for windows
this_dir = os.path.abspath(os.path.dirname(__file__))
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

def getKillerDictionary(killedDict): #makeReciprocalDictionary
    '''
    Take a dictionary of interactions where the first key is the sensitive strain, the second key is the killing strain, and the value is
            the observation frequency
    Return a dictionary of interactions where the first key is the killing strain, the second key is the sensitive strain, and the value is
            the observation frequency
    '''
    killed = killedDict  
    killers = {}
    for strain in killed: 
        for killer in killed[strain]:
            try:
                killers[killer][strain] = killed[strain][killer]
            except:
                killers[killer] = {}
                killers[killer][strain] = killed[strain][killer]

    return killers
     
def combineKillersKilled(killers, killed): #combineNestedDictionaries
    '''
    Join two nested dictionaries on the first key
    Return one nested dictionary
    '''
    traits = {}
    
    def markKillersInKilled(killedDictionary): #markNestedSecondKey
        '''
        Mark the second key in a nested dictionary with an asterisk
        Return a marked nested dictionary
        '''
        killedAst = {}
        for key in killedDictionary:
            killedAst[key] = {}
            for subkey in killedDictionary[key]:
                subkeyAst = subkey + '*'
                killedAst[key][subkeyAst] = killedDictionary[key][subkey]
        return killedAst

    killedAst = markKillersInKilled(killed)
           
    for key in killers:
        traits[key] = dict(killers[key].items()+killedAst[key].items()) 

    return traits
      
def dictionaryToFile(outputFilename, dictionary):
    '''
    '''    
    report = dictionary #trait dictionary     

    lines = []    
    for susceptible in report:  #isolate in report
        for killer in report[susceptible]: #interactor in report[isolate]...
            line = susceptible + "\t" + killer + "\t" + str(report[susceptible][killer]) + "\n"
            lines.append(line)
        
    with open(outputFilename, 'w') as outputFile: #use linesToFile function below
        outputFile.writelines(lines)

def getDistribution(dictionary): #count the hits by superkey of interactions described in the dictionary
    values = []
    for isolate in dictionary:
        count = 0
        for interactor in dictionary[isolate]:
            if dictionary[isolate][interactor] == True:
                count += 1    
        values.append(str(isolate) + ',' + str(count) + "\n")
    return values

def linesToFile(outputFilename, lines):
    with open(outputFilename, 'w') as o:
        o.writelines(lines)

def getDistanceBool(traits):
    '''
    Returns a boolean distance dictionary for isolates in 
    '''
    distanceDict = {}
    isolates = traits.keys()

    for isolate in isolates:    
        distanceDict[isolate] = {}
        for isolate2 in isolates:
            difference = 0
            for trait in traits[isolate]:
                b1 = traits[isolate][trait]
                b2 = traits[isolate2][trait]
                if (b1 == None) or (b2 == None):
                    difference += 0
                elif b1 != b2:
                    difference += 1

            distanceDict[isolate][isolate2] = difference
    
    return distanceDict

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

def createGraph(dictionary):
    nodes = dictionary.keys() #isolates are nodes
    nodePointers = dictionary
    edges = {}
    for eachNode in nodes:
        edges[eachNode] = []
        for killed in nodePointers[eachNode]:
            if nodePointers[eachNode][killed] == True:
                edges[eachNode].append(killed)
    return edges

def identifyLoops(graph):
        
    def identifyUnaryAntagonism(graph):
        uns = []
        for node in graph:
            if node in graph[node]:
                uns.append(node)
        return uns
    
    def identifyBinaryAntagonism(graph):
        bins = []
        for node in graph:
            for node2 in graph[node]:
                if node != node2:
                    if node in graph[node2]:
                        bins.append((node, node2))
        return bins
    
    def identifyTernaryAntagonism(graph):
        terns = []
        for node in graph:
            for node2 in graph[node]:
                if node != node2:
                    for node3 in graph[node2]:
                        if node3 != node2 and node3 != node:
                            if node in graph[node3]:
                                terns.append((node, node2, node3))
        return terns
    
    def identifyQuaternaryAntagonism(graph):
        quats = []
        for node in graph:
            for node2 in graph[node]:
                if node != node2:
                    for node3 in graph[node2]:
                        if node3 != node2 and node3 != node:
                            for node4 in graph[node3]:
                                if node4 != node3 and  node4 != node2 and node4 != node:
                                    if node in graph[node4]:
                                        quats.append((node, node2, node3, node4))
        return quats
    
    return ((identifyUnaryAntagonism(graph)), (identifyBinaryAntagonism(graph)), (identifyTernaryAntagonism(graph)), (identifyQuaternaryAntagonism(graph)))
     
if __name__ == '__main__': #NOT THOROUGHLY TESTED, should be implemented as a wrapper function
    inputFile = 'boolAntagonismAcinetobacter.txt'
    killedDict = loadPairwiseAntagonism(inputFile)
    killerDict = getKillerDictionary(killedDict)
    traits = combineKillersKilled(killedDict, killerDict)
    #dictionaryToFile('traitsAcinetobacter.txt', traits)
    deathCount = getDistribution(killedDict)
    #linesToFile('susceptDistribAcinetobacter.csv', deathCount)
    killCount = getDistribution(killerDict)
    #linesToFile('killDistribAcinetobacter.csv', killCount)
    distance = getDistanceBool(traits)
    dictionaryToMatrix('distanceAcinetobacter.txt', distance)
    dictionaryToMatrix('interactionsAcinetobacter.txt', killedDict)
    graph = createGraph(killerDict)
    abLoops = identifyLoops(graph)
        
    inputFile = 'boolAntagonismVibrio.txt'
    killedDict = loadPairwiseAntagonism(inputFile)
    killerDict = getKillerDictionary(killedDict)
    traits = combineKillersKilled(killedDict, killerDict)
    #dictionaryToFile('traitsVibrio.txt', traits)
    deathCount = getDistribution(killedDict)
    #linesToFile('susceptDistribVibrio.csv', deathCount)
    killCount = getDistribution(killerDict)
    #linesToFile('killDistribVibrio.csv', killCount)
    distance = getDistanceBool(traits)
    dictionaryToMatrix('distanceVibrio.txt', distance)
    dictionaryToMatrix('interactionsVibrio.txt', killedDict)
    graph = createGraph(killerDict)
    vibLoops = identifyLoops(graph)

