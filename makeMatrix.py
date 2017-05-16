from Bio import SeqIO

def getSeqid(file):
    #This function records the names and the sequences within the nexus file.
    # The names of the sequences are placed in the seqid list and will be used to create the first column of the
    # matrix.
    seqid = []
    with open(file, 'r') as fullFile:
        for record in SeqIO.parse(fullFile, "nexus"):
            seqid.append(record.id)
    #Sequence names have been recorded in the list and are returned.
    return seqid
def getSeqs(file):
    #The sequences will be recorded and used to find which of the seqs have the event.
    seqs = []
    with open(file, 'r') as fullFile:
        for record in SeqIO.parse(fullFile, "nexus"):
            seqs.append(record.seq)
    #Sequences have been recorded in the list and are returned.
    return seqs
def getEventNames(file):
    #The gene and the event names will be recorded in the allEvents list. The locations will be recorded in the
    #  eventLocation list. This will be used to find where there is either a dash or a sequence. One or the other
    #  will inform following functions whether a 0 or 1 should be placed in the matrix.
    allEvents = []
    with open(file, 'r') as handle:
        handle = handle.read()
        for lines in handle.splitlines():
            if 'CharSet' in lines:
                words = lines.split(' = ')
                for word in words:
                    if 'CharSet' in word:
                        events = word.split('CharSet ')
                        for event in events:
                            if event == '':
                                pass
                            elif event != '':
                                #All the names of genes and events are recorded.
                                allEvents.append(event)
    return allEvents
def getLocations(file):
    eventNames = getEventNames(file)

def makeMatrix(file):
    seqNames = getSeqid(file)
    seqs = getSeqs(file)
    events = getEventNames(file)
    print(seqNames)
    print(seqs)
    print(events)

makeMatrix('/home/god/Documents/KelcherResearch/testData/LSC1_completed.nex')