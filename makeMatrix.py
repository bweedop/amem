from Bio import SeqIO
import numpy as np

"""The following parser is for NEXUS files exported from the PHYdE program (Muller et al. 2005).
    The events must be labelled as: 'inversion', 'autoinsertion', 'autodeletion', 'autorepeatinsertion',
    'autorepeatdeletion','sharedinsertion', 'shareddeletion', 'sharedrepeatinsertion', 'sharedrepeatdeletion',
    'ambiguous','error', 'review', 'unknown'. Any other events will not be recognized as of now"""


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
    # The gene and the event names will be recorded in the all_events list.
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


def getGenes (file):
    events = getEventNames(file)
    microEvents = ['inversion', 'autoinsertion', 'autodeletion', 'autorepeatinsertion', 'autorepeatdeletion',
                   'sharedinsertion', 'shareddeletion', 'sharedrepeatinsertion', 'sharedrepeatdeletion', 'ambiguous',
                   'error', 'review', 'unknown']
    geneIdentifiers = ['GeneDownstream', 'GeneUpstream']
    pre_genes = [x for x in events if x not in microEvents]
    genes = [y for y in pre_genes if y not in geneIdentifiers]
    return genes


def getLocations(file):
    #Function takes the nexus file and the list of event names as arguments. The locations will be recorded in the
    #  allLocations list.
    microEvents = ['inversion', 'autoinsertion', 'autodeletion', 'autorepeatinsertion', 'autorepeatdeletion',
                   'sharedinsertion', 'shareddeletion', 'sharedrepeatinsertion', 'sharedrepeatdeletion', 'ambiguous',
                   'error', 'review', 'unknown']
    allLocations = []
    with open(file, 'r') as handle:
        handle = handle.read()
        for lines in handle.splitlines():
            if 'CharSet' in lines:
                eventsPlusLocations = lines.split('CharSet ')
                for i in eventsPlusLocations:
                    for a in microEvents:
                        if a in i:
                            locations = i.split(a + ' = ')
                            for j in locations:
                                oneEvent = []
                                eventLocation = j.split(' ')
                                for k in eventLocation:
                                    firstPositions = k.split()
                                    for l in firstPositions:
                                        l = l.replace(';', '')
                                        oneEvent.append(l)
                            allLocations.append(oneEvent)
    return allLocations


def getColumns(file):
    #Allows for the proportional values of each mutational event to be calculated.
    microEvents = ['inversion', 'autoinsertion', 'autodeletion', 'autorepeatinsertion', 'autorepeatdeletion',
                   'sharedinsertion', 'shareddeletion', 'sharedrepeatinsertion', 'sharedrepeatdeletion', 'ambiguous',
                   'error', 'review', 'unknown']
    sorted_events = [x for x in getEventNames(file) if x in microEvents]
    event_values = []
    locations = getLocations(file)
    seq = getSeqs(file)
    list_of_columns = []
    for event_list in locations:
        number_of_columns = 0
        for event in event_list:
            if '-' in event:
                positions = event.split('-')
                number_of_columns += len(range(int(positions[0]), int(positions[1])))
            elif '-' not in event:
                number_of_columns += 1
        list_of_columns.append(number_of_columns)
    while sorted_events:
        for j in list_of_columns:
            percentage = round(j/(len(seq[0])), 6)
            event_values.append(sorted_events.pop(0)+':'+str(percentage))
    return event_values

def dataWrapper (file):
    #Wraps the previous functions in a way that the data is written into a .txt file.
    microEvents = ['inversion', 'autoinsertion', 'autodeletion', 'autorepeatinsertion', 'autorepeatdeletion',
                   'sharedinsertion', 'shareddeletion', 'sharedrepeatinsertion', 'sharedrepeatdeletion', 'ambiguous',
                   'error', 'review', 'unknown']
    sorted_events = [x for x in getEventNames(file) if x in microEvents]
    export = open('data.txt', 'w')
    export.write("Character Sets included in the NEXUS file:\n"+str(getEventNames(file))+'\n')
    current_column = getColumns(file)
    export.write("\nProportional amount of columns pertaining to each mutation type:\n")
    for i in current_column:
        export.write(i+'\n')
    export.write("\nThere are "+str(len(getGenes(file)))+" genes within the file\n. Each gene is listed below:\n"+
                 str(getGenes(file))+'\n')
    export.write("\nThere are "+str(len(getLocations(file)[sorted_events.index('inversion')]))+" inversions.\n"+
                 "Positions of inversions:\n"+str(getLocations(file)[sorted_events.index('inversion')])+'\n')
    export.write("\nThere are "+str(len(getLocations(file)[sorted_events.index('ambiguous')]))+
                 " ambiguous regions.\n"+"Positions of ambiguous regions:\n"+
                 str(getLocations(file)[sorted_events.index('ambiguous')])+'\n')
    export.write("\nThere are "+str(len(getLocations(file)[sorted_events.index('review')]))+
                 " regions needing to be reviewed.\n"+ "Positions needing to be reviewed:\n"+
                 str(getLocations(file)[sorted_events.index('review')])+'\n')
    export.write("\nThere are "+str(len(getLocations(file)[sorted_events.index('unknown')]))+ " unknown events.\n"+
                 "Positions of unknown events:\n"+str(getLocations(file)[sorted_events.index('unknown')])+'\n')
    export.close()