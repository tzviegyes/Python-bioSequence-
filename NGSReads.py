import numpy as np
from typing import List

class NgsReads:
    """A class of NGS DNA Reads"""
    def __init__(self, name:str, sequence:str, quality_sequence:str=''):
        """constructor for NGS Read object.
        :param name: name of reads as apears in fq file
        :param sequence: the DNA sequence of letters ACTG
        :param quality_sequence: quality sequence of Phred33"""
        self.name = name
        self.sequence = sequence
        self.string_quality = quality_sequence
        self.quality_sequence = np.array([ord(letter)-33 for letter in quality_sequence])
        self.low_quality_kmers = 0

    def __all_kmers(self, K:int, mkq:float):
        """recives read and K and minimum kmer quality returns all Kmers of length K that their quality is high enough"""
        all_kmers = []
        length = len(self.sequence)
        start_position = 0
        if len(self.quality_sequence) > 0 and len(self.quality_sequence) == length:
            while start_position <= (length - K):
                kmer_quality = np.average(self.quality_sequence[start_position: start_position + K])
                if kmer_quality >= mkq:
                    all_kmers.append(self.sequence[start_position: start_position + K])
                else:
                    self.low_quality_kmers += 1
                start_position += 1
        else:
            while start_position <= (length - K):
                all_kmers.append(self.sequence[start_position: start_position + K])
                start_position += 1
        return all_kmers

    def __read_quality_maen(self) -> float:
        """returns the quality of a given sequence"""
        if len(self.quality_sequence) > 0:
            quality = np.average(self.quality_sequence)
            return quality
        else:
            return 256

    def __reverse_complement(self):
        """returns the reversee complement of a read"""
        new_sequence = ''
        for letter in self.sequence[::-1]:
            if letter == 'A':
                new_sequence += 'T'
            if letter == 'T':
                new_sequence += 'A'
            if letter == 'G':
                new_sequence += 'C'
            if letter == 'C':
                new_sequence += 'G'
        reversed_read = NgsReads(self.name + 'r', new_sequence, self.string_quality[::-1])
        return reversed_read

    def get_name(self) -> str:
        """geter function for read's name"""
        return self.name

    def get_sequence(self) -> str:
        """getter function for the read's DNA sequence"""
        return self.sequence

    def get_read_qulity(self) -> float:
        """getter function for the read's quality mean"""
        return self.__read_quality_maen()

    def get_all_kmers(self, K:int, mkq:float=0) -> List[str]:
        """getter function for all K-mers in read"""
        return self.__all_kmers(K, mkq)

    def get_low_qulity_count(self) -> int:
         """getter function for the counter of K-mers that are skipped dou to low quality"""
         return self.low_quality_kmers

    def get_reverse_complement(self):
        return self.__reverse_complement()

