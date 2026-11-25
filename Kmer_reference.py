from typing import Dict, List, Union, Any, Optional
from Reference import Reference


class KmerReference:
    """The K-mer reference is the main database the psuedo alginment uses.
    It is implemented as a python dictionary (hash-table) of dictionaries.
    key: K-mer, value: {key: genome, value: list of positions the K-mer apears in the genome}
    For example: {"AGT": {dog: [2, 8, 11], cat: [0, 10]}}
    meaning the sequence "AGT" apears three times in the genome of a dog and twice in the genome of a cat.
    """
    def __init__(self, data:Any, K=0, filter:Union[bool, float] = False):
        """
        A constructor for a KmerReference object.
        :param data: could get a list of References object - the genomes the database will include.
                    or a K-mer reference summary.
        :param K: lenght of Kmers
        """
        self._KmerDB: Dict = {}
        self.K = K
        self.filtering_dict = {}
        self.already_filtered = False
        self.original_gen_order:List[Union[str,Reference]] = []
        self.filter = filter

        if isinstance(data, list):
            self.lenght = len(data)
            self.__genomes = data
            self.true_if_list = True

        if isinstance(data, dict):
            self.true_if_list = False
            self.data = data["reference_database"]
            self.filtering_dict = data["filter_stats"]
            self.__genomes = []

    def __valid_input(self) -> bool:
        """checks all genomes are valid, all have the same K"""
        if self.true_if_list:
            for genome in self.__genomes:
                if not genome.get_is_input_valid() or genome.get_K() != self.K:
                    return False
            return True
        else:
            return True

    def __construct_from_list(self) -> Dict[str, Dict[Reference, List[int]]]:
        """constructs the database from list of References"""
        if self._KmerDB == {}:
            for genome in self.__genomes:
                kmer_dict_genome = genome.get_reference_kmer_dict()
                for Kmer in kmer_dict_genome:
                    if Kmer in self._KmerDB:
                        self._KmerDB[Kmer].update({genome: kmer_dict_genome[Kmer]})
                    else:
                        self._KmerDB.update({Kmer: {genome: kmer_dict_genome[Kmer]}})
            if self.original_gen_order == []:
                self.original_gen_order = [str(gen) for gen in self.__genomes]

            if self.filter and self.already_filtered == False:
                self.already_filtered = True
                self.__greedy_filtering()
                self._KmerDB = {}
                self.__construct_from_list()
        return self._KmerDB

    def __list_of_genomes(self) -> List[Union[Reference, str]]:
        """re-constructs list of genomes from a pre-built database"""
        if not self.true_if_list:
            list_of_gens = []
            for sub_dict in self.data.values():
                for gen in sub_dict.keys():
                    if gen not in list_of_gens:
                        list_of_gens.append(gen)
            self.__genomes = list_of_gens
        return self.__genomes

    def __K_from_dict(self) -> int:
        """returns the K-mer lenght, wheather given from user or infered from database"""
        if not self.true_if_list:
            self.K = len(list(self.data.keys())[0])
        return self.K

    @property
    def KmerDB(self) -> Dict[str, Dict[Union[Reference,str], List[int]]]:
        """calls construct_from_list if data is list of references,
        otherwise assigns the DB atriburte to be the input data"""
        if self._KmerDB != {}:
            return self._KmerDB
        if not self.__valid_input():
            print("invalid K-mer input")
            exit()
        if self.true_if_list:
            return self.__construct_from_list()
        else:
            self._KmerDB = self.data
            return self.data

    def __summerize_genomes(self) -> Dict[Union[str,Reference], Dict[str, int]]:
        """produces a summary of all genomes in database"""
        database = self.KmerDB
        genomes_dict = {genome: {"total_bases": self.__K_from_dict() - 1,
                                 "unique_kmers": 0,
                                 "multi_mapping_kmers": 0} for genome in self.__list_of_genomes()}
        for kmer in database:
            sub_dict = database[kmer]
            if len(sub_dict) == 1: # unique
                genomes_dict[list(sub_dict.keys())[0]]["unique_kmers"] += 1
                genomes_dict[list(sub_dict.keys())[0]]["total_bases"] += len(list(sub_dict.values())[0])
            else:
                for gen,place_list in sub_dict.items():
                    genomes_dict[gen]["total_bases"] += len(place_list)
                    genomes_dict[gen]["multi_mapping_kmers"] += 1
        return genomes_dict

    def __dict_for_filter(self) -> Dict[Reference, Dict[str, int]]:
        """construct a dictionary with all genomes and the information needed for filtering.
        Then sorts by specific count, total count and lenght"""
        genome_summary = self.__summerize_genomes()
        for gen,gen_summary in genome_summary.items():
            self.filtering_dict.update({gen: {"unique_kmers": gen_summary["unique_kmers"],
                                                   "total_kmers": gen_summary["unique_kmers"] + gen_summary["multi_mapping_kmers"],
                                                   "genome_length": gen_summary["total_bases"]}})
        sorted_dictionary = dict(sorted(self.filtering_dict.items(),
                                        key=lambda item: (item[1]["unique_kmers"], item[1]["total_kmers"],
                                                          item[1]["genome_length"])))
        self.filtering_dict = sorted_dictionary
        return self.filtering_dict


    def __greedy_filtering(self) -> List[Reference]:
        """filters over the genoms and removes similar genomes.
        keeps track of statistics"""
        genome_list = list(self.__dict_for_filter().keys())
        out_index = 0
        while out_index < len(genome_list):
            gen1 = genome_list[out_index]
            gens1_kmers = gen1.get_reference_kmer_dict()
            if isinstance(gens1_kmers, bool):
                out_index += 1
                continue
            else:
                kmers1 = set(gens1_kmers.keys())
            in_index = out_index + 1
            while in_index < len(genome_list):
                gen2 = genome_list[in_index]
                gen2_kmers = gen2.get_reference_kmer_dict()
                if isinstance(gen2_kmers, bool):
                    in_index += 1
                    continue
                else:
                    kmers2 = set(gen2_kmers.keys())
                similar_score = (len(kmers1 & kmers2) / min(len(kmers1), len(kmers2)))
                if similar_score > self.filter:
                    self.filtering_dict[gen1] = ({"kept": "no"} | (self.filtering_dict[gen1])
                                                 | {"similar_to": str(gen2), "similarity_score": similar_score})
                    self.__genomes.remove(gen1)
                    # genome_list.remove(gen1)
                    break
                in_index += 1
            else:
                self.filtering_dict[gen1] = ({"kept": "yes"} | self.filtering_dict[gen1]
                                             | {"similar_to": 'NA', "similarity_score": 'NA'})
            out_index += 1
        return self.__genomes

    def get_genomes(self) -> List[Union[str, Reference]]:
        """geter function for list of genomes """
        genomes = self.__list_of_genomes()
        return genomes

    def get_K(self) -> int:
        """getter function for the K-mer length"""
        k = self.__K_from_dict()
        return k

    def get_kmer_referenceDB(self) -> Dict[Union[str,Reference], Dict[str, List[int]]]:
        """getter function for the database"""
        data_base = self.KmerDB
        return data_base

    def get_valid_input(self) -> bool:
        """getter function to knoe if the input is valid"""
        return self.__valid_input()

    def get_filter_stats(self):
        if not self.filter:
            return self.filtering_dict
        else:
            database = self.KmerDB
            stats_for_printing = {str(gen): stats for gen,stats in self.filtering_dict.items()}
            sorted_by_original_order = {gen: stats_for_printing[gen] for gen in self.original_gen_order}
            return sorted_by_original_order

    def get_database_for_printing(self) -> Dict[str, Dict[str, List[int]]]:
        """getter function for the database in a printable form"""
        database = self.KmerDB
        dict_for_print = {}
        for kmer in database:
            temp_dict = {}
            for gen,places in database[kmer].items():
                temp_dict.update({str(gen): places})
            dict_for_print.update({kmer: temp_dict})
        return dict_for_print

    def get_genome_summery(self) -> Dict[str, Dict[str, int]]:
        """getter function for the genome summary"""
        gen_summary = {str(gen): val for gen, val in self.__summerize_genomes().items()}
        return gen_summary
