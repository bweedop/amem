import xlrd
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_alphabet
from Bio.SeqRecord import SeqRecord

def loadSheet(file):
    #loads the .xlsx file into memory using the xlrd module.
    char_data = xlrd.open_workbook(file)
    sheet = char_data.sheet_by_index(0)

    return sheet

def extractLocations(file):
    #Extracts the event locations from the 'columns' column in the .xlsx file.
    char_data = loadSheet(file)
    a_z = list(map(chr, range(97, 123)))
    a_z.extend(list(map(chr, range(65, 91))))
    cells = [[str(char_data.cell_value(c, 1))] for c in range(char_data.nrows)]
    for i in cells:
        if str(i)[2] in a_z:
            cells.remove(i)
    cells.pop(0)

    return cells

def extractType(file):
    #Extracts the type of each event from the 'type' column in the .xlsx file.
    char_data = loadSheet(file)
    cells = [[str(char_data.cell_value(c, 3))] for c in range(char_data.nrows)]
    for i in cells:
        if i == ['']:
            cells.remove(i)
    cells.pop(0)

    return cells


def extractTaxa(file):
    #Extracts the taxa names from the 'taxa' column in the .xlsx file.
    char_data = loadSheet(file)
    cells = [[str(char_data.cell_value(c, 8))] for c in range(char_data.nrows)]
    for i in cells:
        if i == ['']:
            cells.remove(i)
    cells.pop(0)
    only_taxa = []
    for j in cells:
        for k in range(len(j)):
            a = [x.split() for x in j[k].replace(')', '').split(' (?')]
            only_taxa.append(a)
    return only_taxa


def makeMatrix(file, file_path):
    """Wrapper function that creates the character matrix and writes the matrix to a .csv file. The matrix is
    initialized as a numpy array filled with zeros. A one is inserted to the array where events have been determined
    to be present. A '?' is inserted for taxa where an event may have occurred but cannot be verified. The array is
    converted to a pandas dataframe, given row names (sequence names) and column names (1- #of events), and
    written to a .csv file. The user must give two arguments: (1) the .xlsx file data is being pulled from and
    (2) the path and file name desired for output."""
    locations = extractLocations(file)
    taxa_lists = extractTaxa(file)
    seqs = []
    for i in taxa_lists:
        for j in i:
            for k in j:
                if k not in seqs:
                    seqs.append(k)
    char_matrix = np.zeros((len(seqs), len(locations)), dtype=object)
    event_position = 0
    for a in taxa_lists:
        event_position += 1
        if len(a) == 1:
            for b in a:
                for c in b:
                    char_matrix[seqs.index(c), event_position-1] = 1
        elif len(a) != 1:
            for b in a[0]:
                char_matrix[seqs.index(b), event_position - 1] = 1
            for d in a[1]:
                char_matrix[seqs.index(d), event_position - 1] = '?'
    df = pd.DataFrame(char_matrix, index=seqs, columns=[x+1 for x in range(len(locations))])
    df.to_csv(file_path)
