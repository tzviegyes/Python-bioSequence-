# from xml.etree.ElementTree import indent

from Reference import Reference
from NGSReads import NgsReads
import csv
import pickle
from typing import List, Any

def make_list_of_genomes(file_name:str, K:int) -> List[Reference]:
    """Recives a .fa file of reference genomes.
    Returns list of Reference objects"""
    if K < 1:
        print("K-mer size needs to be at least 1.\n")
        exit(1)
    try:
        genomes_list = []
        with open(file_name, 'r') as reference_file:
            reference_reader = csv.reader(reference_file)
            current_header = None
            current_sequence = ''

            for row_list in reference_reader:
                row = row_list[0]
                if row.startswith('>'):
                    # add previous genome
                    if current_header is not None:
                        genomes_list.append(Reference(current_header, current_sequence, K))
                    # start new genome
                    current_header = row[1:]
                    current_sequence = ''
                else:
                    current_sequence += row
            if current_header is not None:
                genomes_list.append(Reference(current_header, current_sequence, K)) # add last genome
        if genomes_list == []:
            print("Genome reference file seems to be empty or not in .fa format\n")
            exit(1)
        return genomes_list
    except (FileExistsError,FileNotFoundError):
        print("Error. File not found")
        exit(1)

def make_list_of_reads(file_name:str, with_quality=False) -> List[NgsReads]:
    """Recives a .fq file of NGS reads.
    Returns list tuples (name, sequence)"""
    try:
        reads_list = []
        with open(file_name, 'r') as reads_file:
            reads_reader = csv.reader(reads_file)
            for index, row_list in enumerate(reads_reader):
                row = row_list[0]
                if index % 4 == 0:
                    current_header = row[1:]
                if index % 4 == 1:
                    current_row = row
                if index % 4 == 2:
                    continue
                if index % 4 == 3:
                    if with_quality is True:
                        reads_list.append(NgsReads(current_header, current_row, row))
                    else:
                        reads_list.append(NgsReads(current_header, current_row))
        if reads_list == []:
            print("NGS reads file seems to be empty or not in .fastq format.\n")
            exit(1)
        return reads_list

    except (FileExistsError, FileNotFoundError):
        print("Error. File not found")
        exit(1)

def write_data_to_file(file_name:str, data:Any) -> bool:
    """
    :param file_name: name of file to write the kbd file
    :return: True if completed tesk
    """
    try:
        with open(file_name, "wb") as file:
            pickle.dump(data, file)
        return True
    except (FileExistsError, FileNotFoundError):
        print("Error. could'nt write to file or path")
        exit(1)

def read_data_file(file_name:str) -> Any:
    """recives a compressed pickle file and returns the data in the file"""
    try:
        with open(file_name, "rb") as file:
            loaded_data = pickle.load(file)
            if not loaded_data:
                print("Problem accured when reading file.\n")
                exit(1)
            return loaded_data
    except (FileExistsError, FileNotFoundError):
        print("Error. File not found\n")
        exit(1)

