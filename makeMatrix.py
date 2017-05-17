from Bio import SeqIO


def getSeqid(file):
    # This function records the names and the sequences within the nexus file.
    # The names of the sequences are placed in the seq_id list and will be used to create the first column of the
    # matrix.
    seq_id = []
    with open(file, 'r') as fullFile:
        for record in SeqIO.parse(fullFile, "nexus"):
            seq_id.append(record.id)
    # Sequence names have been recorded in the list and are returned.
    return seq_id


def getSeqs(file):
    # The sequences will be recorded and used to find which of the seqs have the event.
    seqs = []
    with open(file, 'r') as fullFile:
        for record in SeqIO.parse(fullFile, "nexus"):
            seqs.append(record.seq)
    # Sequences have been recorded in the list and are returned.
    return seqs


def getEventNames(file):
    # The gene and the event names will be recorded in the all_events list.
    all_events = []
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
                                # All the names of genes and events are recorded.
                                all_events.append(event)
    return all_events


def getLocations(file, event):
    # Function takes the nexus file and the event name as arguments. The locations will be recorded in the
    #  all_Location list.
    all_locations = []
    with open(file, 'r') as handle:
        handle = handle.read()
        for lines in handle.splitlines():
            if 'CharSet' in lines:
                events_plus_locations = lines.split('CharSet ')
                for i in events_plus_locations:
                    if event in i:
                        locations = i.split(event + ' = ')
                        for j in locations:
                            event_location = j.split(' ')
                            for k in event_location:
                                first_positions = k.split()
                                for l in first_positions:
                                    l = l.replace(';', '')
                                    all_locations.append(l)
    return all_locations


