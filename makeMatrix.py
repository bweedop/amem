from Bio import SeqIO

def getSeqid(file):
    """Record names and sequences within a nexus file
    
    The names of the sequences are placed in the seqid list and will
    be used to create the first column of the matrix.
    """
    with open(file, 'r') as fullFile:
        output = [record.id for record in SeqIO.parse(fullFile, "nexus")]
    
    return output

def getSeqs(file):
    #The sequences will be recorded and used to find which of the seqs have the event.
    with open(file, 'r') as fullFile:
        seqs = [record.seq for record in SeqIO.parse(fullFile, "nexus")]
    #Sequences have been recorded in the list and are returned.
    return seqs

def getEventNames(file):
    #The gene and the event names will be recorded in the allEvents list.
    allEvents = []
    with open(file, 'r') as handle:
        handle = handle.read()
        for lines in handle.splitlines():
            if 'CharSet' in lines:
                words = lines.split(' = ')
                for word in words:
                    if 'CharSet' in word:
                        allEvents.append(word.replace('CharSet ', ''))
    
    return allEvents

def getLocations(file, events):
    #Function takes the nexus file and the list of event names as arguments. The locations will be recorded in the
    #  eventLocation list.
    allLocations = []
    with open(file, 'r') as handle:
        handle = handle.read()
        for lines in handle.splitlines():
            if 'CharSet'in lines:
                eventsPlusLocations = lines.split('CharSet ')
                currentEvent = events.pop(0)
                for i in eventsPlusLocations:
                    if currentEvent in i:
                        locations = i.split(currentEvent+' = ')
                        for j in locations:
                            eventLocation = j.split(' ')
                            for k in eventLocation:
                                firstPositions = k.split()
                                for l in firstPositions:
                                    allLocations.append(l.replace(';',''))
    return allLocations
