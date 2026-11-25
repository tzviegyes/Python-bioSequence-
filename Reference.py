from typing import Dict, List, Union

class Reference:
    """A Reference is a DNA sequence of a known organizm"""
    ALLOWED_LETTERS = "ACGTN"
    def __init__(self, name:str, sequence:str, K:int) -> None:
        """
        A constructor for a Reference object.
        :param name: A string representing the Reference's name.
        :param sequence: The sequence of DNA. for example "AAGTCA"
        :param K: lenght of Kmers
        """
        self._ref_kmer_dict:Dict = {}
        self.name = name
        self.sequence = sequence
        self.lengh = len(sequence)
        self.K = K

    def __str__(self) -> str:
        return str(self.name)

    def __valid_input_name(self) -> bool:
        """returns True if name input is string"""
        return isinstance(self.name, str)

    def __valid_input_sequence(self) -> bool:
        """returns True if sequence input is string, and contains only AGTCN letters"""
        if not isinstance(self.sequence, str):
            return False
        for nucleotid in self.sequence:
            if nucleotid not in self.ALLOWED_LETTERS:
                return False
        return True

    def __valid_input_K(self) -> bool:
        """returns True if K input is positive interger and shorter than sequence"""
        if not isinstance(self.K, int):
            return False
        return self.lengh >= self.K

    def _valid_all_input(self) -> bool:
        """checks all input is valid"""
        if not self.__valid_input_name() or not self.__valid_input_sequence():
            return False
        return self.__valid_input_K()

    @property
    def ref_kmer_dict(self) -> Union[Dict[str, List[int]], bool]:
        """builds a dictionary of Kmer: list of positions the Kmer apears in the sequence.
        contains only Kmers that do not contain N"""
        if not self._valid_all_input():
            return False
        if self._ref_kmer_dict == {}: # Check it wasn't built yet
            start_position = 0 if 'N' not in self.sequence[:self.K] else self.sequence[:self.K].index('N') + 1
            while start_position <= (self.lengh - self.K):
                if self.sequence[start_position + self.K - 1] == 'N':
                    start_position += self.K
                    continue
                if self.sequence[start_position] == 'N':
                    start_position += 1
                    continue
                temp_sequence = self.sequence[start_position: start_position + self.K]
                if self._ref_kmer_dict is not None and temp_sequence in self._ref_kmer_dict:
                    self._ref_kmer_dict[temp_sequence].append(start_position)
                else:
                    self._ref_kmer_dict[temp_sequence] = [start_position]
                start_position += 1
        return self._ref_kmer_dict

    def get_K(self) -> int:
        """getter function for K-mer length"""
        return self.K

    def get_sequence(self) -> str:
        """getter function for the genome's DNA sequence"""
        return self.sequence

    def get_is_input_valid(self) -> bool:
        """returns True if the input is valid"""
        return self._valid_all_input()

    def get_reference_kmer_dict(self) -> Union[Dict[str, List[int]], bool]:
        """getter function for the K-mer positions dictionary"""
        dicitionary = self.ref_kmer_dict
        return dicitionary

