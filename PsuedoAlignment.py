from Kmer_reference import KmerReference
from Reference import Reference
from NGSReads import NgsReads
from typing import List, Dict, Tuple, Union, Optional

class PsuedoAlginment:
    """The Psuedo alginment algorithm"""
    ALLOWED_LETTERS = "ACGTN"
    def __init__(self, ngs_reads:List[NgsReads], kmer_reference:KmerReference,  m=1, p=1, min_read_quality:float=256,
                 min_kmer_quality:float=0, max_gen=-1, check_reverse = False):
        """
        A constructor for a PsuedoAlginment object.
        :param ngs_reads: the DNA samples the algorithm will map.
        :param kmer_reference: the database of known genomes the algorithm is looking in.
        :param m: minimal difference between the highest and second uniquely mappaed genomes.
        :param p: maximal diference between total count of initaly uniquely mapped to other genomes.
        """
        self.ngs_current_read = ngs_reads[0]
        self.reads = ngs_reads
        self.kmer_referenceDB = kmer_reference.get_kmer_referenceDB()
        self.kmer_reference_genomes = kmer_reference.get_genomes()
        self._genome_stats:Dict = {}
        self._mapped_genome:Union[str,Reference] = ""
        self.K = kmer_reference.get_K()
        self.p = p
        self.m = m
        self.min_read_quality = min_read_quality
        self.min_kmer_quality = min_kmer_quality
        self.max_genomes = max_gen
        self.__skipped_reads = 0
        self.__skipped_low_quality_kmers = 0
        self.__skipped_high_mapp_kmers = 0
        self.temp_high_mapp_counter = 0
        self.temp_low_qulity_kmers = 0
        self.reads_stats = {"unique_mapped_reads":0 ,"ambiguous_mapped_reads":0, "unmapped_reads":0}
        self.kmers_mapping = ([''],[''])
        self.to_check_reverese = check_reverse

    def __valid_input(self) -> bool:
        """checks: K shorter than sequence, sequence are letters ACTGN.
        p and m are integers
        assumes it gets a valid database"""
        if not isinstance(self.ngs_current_read, NgsReads) or not isinstance(self.m, int) or not isinstance(self.p, int) or not isinstance(self.K, int):
            return False
        for letter in self.ngs_current_read.get_sequence():
            if letter not in self.ALLOWED_LETTERS:
                return False
        if self.K > len(self.ngs_current_read.get_sequence()):
            return False
        return True

    @property
    def build_ambiguous_gen_map(self) -> Dict[str, Dict[str, int]]:
        """initiates a mapping of each genome in the references to 0"""
        if self._genome_stats == {}:
            self._genome_stats = {str(genome): {"unique_reads":0, "ambiguous_reads":0} for genome in self.kmer_reference_genomes}
            return self._genome_stats
        else:
            return self._genome_stats

    def __are_kmers_specific(self) -> Tuple[List[str], List[str]]:
        """returns a list all specific Kmers and a list of all unspecific Kmers.
        doesn't include kmers that do not map to any genome."""
        all_kmers = self.ngs_current_read.get_all_kmers(self.K, self.min_kmer_quality)
        self.temp_low_qulity_kmers = self.ngs_current_read.get_low_qulity_count()
        specific_kmers = []
        un_specific_kmers = []
        self.temp_high_mapp_counter = 0
        for kmer in all_kmers:
            if kmer in self.kmer_referenceDB:
                if len(self.kmer_referenceDB[kmer]) > 1:
                    if len(self.kmer_referenceDB[kmer]) > self.max_genomes > -1:
                        self.temp_high_mapp_counter += 1
                    else:
                        un_specific_kmers.append(kmer)
                else:
                    specific_kmers.append(kmer)
        self.kmers_mapping = (specific_kmers, un_specific_kmers)
        return self.kmers_mapping

    def __is_uniquly_mapped(self) -> Union[bool, Union[Reference, str]]:
        """returns True if initialy the read is uniquly mapped"""
        specific_kmers = self.kmers_mapping[0]
        database = self.kmer_referenceDB
        first_genome = list(database[specific_kmers[0]].keys())[0]
        for kmer in specific_kmers:
            if database[kmer].keys != first_genome:
                return False
        self._mapped_genome = first_genome
        return first_genome

    def __genome_any_kmer_count(self, kmers:List[str]) -> Dict[Union[Reference, str], int]:
        """returns a dictionary of all genomes from database and the number of specific or unspecific kmers linked to genome."""
        genom_dict:Dict[Reference | str, int] = {genome: 0 for genome in self.kmer_reference_genomes}
        database = self.kmer_referenceDB
        for kmer in kmers:
            temp_dict = database[kmer]
            for gen,positions in temp_dict.items():
                genom_dict[gen] += len(positions)
        return genom_dict

    def __genome_specific_kmer_count(self) -> Dict[Union[Reference, str], int]:
        """returns a dictionary of all genomes from database and the number of specific kmers linked to genome."""
        specific_kmers = self.kmers_mapping[0]
        return self.__genome_any_kmer_count(specific_kmers)

    def __find_max_and_second(self) -> List[Tuple[Union[Reference, str],int]]:
        """finds the genome with the most unique Kmers mappad to it and the second most"""
        genom_dict = self.__genome_specific_kmer_count()
        max_genome = max(genom_dict, key=lambda gen: genom_dict[gen])
        max_value = genom_dict[max_genome]
        genom_dict.pop(max_genome)
        if len(genom_dict) == 1:
            return [(max_genome, max_value)]

        second_genome = max(genom_dict, key=lambda gen: genom_dict[gen])
        second_value = genom_dict[second_genome]
        genom_dict.update({max_genome: max_value})
        return [(max_genome, max_value),(second_genome, second_value)]

    def __genome_un_specific_kmer_count(self) -> Dict[Union[Reference, str],int]:
        """returns a dictionary of all genomes from database and the number of unspecific kmers linked to genome."""
        un_specific_kmers = self.kmers_mapping[1]
        return self.__genome_any_kmer_count(un_specific_kmers)

    def __is_highest_enough(self) -> Union[Union[Reference, str], bool]:
        """Find the difference between the counts of the highest and second-highest unique k-mer count.
        Returns the genome if the difference is greater or equal to m, if not returns False"""
        two_max = self.__find_max_and_second()
        if len(two_max) == 1:
            self._mapped_genome = two_max[0][0]
            return two_max[0][0]
        if (two_max[0][1] - two_max[1][1]) >= self.m:
            self._mapped_genome = two_max[0][0]
            return two_max[0][0]
        else:
            return False

    def __incriment_ambiguous(self, genomes:Dict[Union[Reference, str], int], threshold=1) -> List[str]:
        """recives a dictionary of {genomes:integer}
        param threshold: determines from what count to incriment genome. for specific the threshold is 0 (2.2.2.6).
        for total count the threshold is MapCount (2.3.2.4).
        returns a list of all genoms ambiguousely mapped"""
        mapped_genoms = []
        subset_genomes = dict(filter(lambda x: x[1] >= threshold, genomes.items()))
        for gen in subset_genomes:
            gen = str(gen)
            self._genome_stats[gen]["ambiguous_reads"] += 1
            mapped_genoms.append(gen)
        return mapped_genoms

    def __total_kmers_for_genomes(self) -> Dict[Union[Reference, str], int]:
        """recives all genomes and their count of specific and unspecific kmers from read.
        returns a dictionay of genomes and total Kmers mapped to genome"""
        specific_gens = self.__genome_specific_kmer_count()
        un_specific_gens = self.__genome_un_specific_kmer_count()
        total_dict = un_specific_gens.copy()  # dict(filter(lambda x: x[1] != 0, specific_gens.items()))
        for gen in total_dict:
            if gen in specific_gens:
                total_dict[gen] += specific_gens[gen]
        return total_dict

    def __validate_uniquley_mapped(self) -> bool:
        """validats that the initial assignmen of read as uniquley mappad is correct
        according to the unspecific Kmers too. returns True if validated.
        Check if MaxCount - MapCount > p. if True: change mapping to ambiguous and incriment count and return False"""
        total_kmers = self.__total_kmers_for_genomes()
        max_count = total_kmers[max(total_kmers, key=lambda kmer: total_kmers[kmer])]
        if self._mapped_genome != "":
            map_count = total_kmers[self._mapped_genome]
        else:
            return False
        if max_count - map_count > self.p:
            return False
        return True

    def __choose_read_orientation(self) -> Tuple[List[str],List[str]]:
        """if user asked to check reverse complement, this function checks the specific K-mer count for both directions
         and chooses the orientation with higher count.
         anyway, returns the specific and not specific K-mers of the read."""
        if self.to_check_reverese:
            forward_count = self.__are_kmers_specific()
            forward_high_kmer = self.temp_high_mapp_counter
            forward_low_quality_kmers = self.temp_low_qulity_kmers
            forward_read = self.ngs_current_read
            self.ngs_current_read = self.ngs_current_read.get_reverse_complement()
            reverse_count = self.__are_kmers_specific()
            reverse_high_kmer = self.temp_high_mapp_counter
            reverse_low_quality_kmers = self.temp_low_qulity_kmers
            if len(reverse_count[0]) > len(forward_count[0]):
                self.__skipped_low_quality_kmers += forward_low_quality_kmers
                self.__skipped_high_mapp_kmers += reverse_high_kmer
                self.kmers_mapping = reverse_count
                return reverse_count
            else:
                self.__skipped_low_quality_kmers += reverse_low_quality_kmers
                self.ngs_current_read = forward_read
                self.__skipped_high_mapp_kmers += forward_high_kmer
                self.kmers_mapping = forward_count
                return forward_count
        else:
            specific_and_not = self.__are_kmers_specific()
            self.__skipped_low_quality_kmers += self.temp_low_qulity_kmers
            self.__skipped_high_mapp_kmers += self.temp_high_mapp_counter
            self.kmers_mapping = specific_and_not
            return specific_and_not

    def __psuedo_align_read(self) -> Optional[Tuple[bool,List[Union[str,Reference]]]]:
        """runs the algorithm using the functions defined above
        returns False if read is unmapped, and an empty list of genomes it is mapped to.
        returns True if read is uniqly mapped and a list of the only genome the read is mapped to.
        returns False if read is ambiguously mapped and a list of genomes it is mapped to."""
        self.build_ambiguous_gen_map
        if not self.__valid_input():
            print("invalid input")
            return None
        # step 1
        specific_kmers, un_specific_kmers = self.__choose_read_orientation()  # 1.2

        # step 2
        if specific_kmers == []: # 2.1
            self.reads_stats["unmapped_reads"] += 1
            return False,[]

        if not self.__is_uniquly_mapped(): # 2.2.2
            if not self.__is_highest_enough(): # 2.2.2.6
                self.reads_stats["ambiguous_mapped_reads"] += 1
                return False, self.__incriment_ambiguous(self.__genome_specific_kmer_count(), 1)

        validated = self.__validate_uniquley_mapped() # 2.3.2.4
        if not validated:
            self.reads_stats["ambiguous_mapped_reads"] += 1
            return False, self.__incriment_ambiguous(self.__total_kmers_for_genomes(), self.__total_kmers_for_genomes()[self._mapped_genome])
        if validated:
            self.reads_stats["unique_mapped_reads"] += 1
            if self._mapped_genome != "":
                self._genome_stats[str(self._mapped_genome)]["unique_reads"] += 1
            return True, [self._mapped_genome]

    def psuedo_algin_reads_list(self) -> Optional[Dict[str, Tuple[bool, List[str]]]]:
        """returns a dictionary of all reads and their mapping"""
        mapping_dict = {}
        for read in self.reads:
            self.kmers_mapping = ([''], [''])
            if read.get_read_qulity() >= self.min_read_quality or self.min_read_quality == 256:
                self.ngs_current_read = read
                mapping_dict.update({read.get_name(): self.__psuedo_align_read()})
            else:
                self.__skipped_reads += 1
        return mapping_dict

    def get_reads_stats(self) -> Dict[str, int]:
        """getter functions for the reads' statistics"""
        self.psuedo_algin_reads_list()
        if self.min_read_quality != 256:
            self.reads_stats.update({"filtered_quality_reads": self.__skipped_reads})
        if self.min_kmer_quality != 0:
            self.reads_stats.update({"filtered_quality_kmers": self.__skipped_low_quality_kmers})
        if self.max_genomes != -1:
            self.reads_stats.update({"filtered_hr_kmers": self.__skipped_high_mapp_kmers})
        return self.reads_stats

    def get_mapping_for_read(self) -> Optional[Tuple[bool,List[Union[str, Reference]]]]:
        """getter function for read's mapping"""
        return self.__psuedo_align_read()

    def get_genome_stats(self) -> Dict[str, Dict[str, int]]:
        """getter function for the genome statistics"""
        return self._genome_stats

    def get_psuedo_algin_reads_list(self) -> Optional[Dict[str, Tuple[bool, List[str]]]]:
        """getter function for dictionary of all reas and their mapping"""
        mapping_dict = self.psuedo_algin_reads_list()
        return mapping_dict

    def get_skipped_low_quality_kmer(self) -> int:
        return self.__skipped_low_quality_kmers

    def get_skiped_reads(self) -> int:
        skipped = self.__skipped_reads
        return skipped

    def get_skipped_kmers_high_mapping(self) -> int:
        return self.__skipped_high_mapp_kmers

class AlignmentData:
    """A class for an alignment summary"""

    def __init__(self, align_summary:Dict[str,Tuple[bool,List[str]]]):
        """initiator for the AlignmentData class. recives a summary:
        dictionary of key: read, value:mapping.
        mapping: (True, genome) if uniquly mapped
                (False []) if unmapped
                (False, [genomes]) if ambiguously mapped"""
        self.summary = align_summary
        self.reads_stats = {"unique_mapped_reads":0 ,"ambiguous_mapped_reads":0,"unmapped_reads":0}
        self.genome_stats:Dict[str, Dict[str, int]] = {}

    def __incriment_genome(self, gen:str, true_if_uniqe:bool) -> None:
        """recives a genome and True if it uniquly mapped, False if ambiguous.
        incriments the counter of genome"""
        if gen not in self.genome_stats:
            self.genome_stats.update({gen: {"unique_reads":0,"ambiguous_reads":0}})
        if true_if_uniqe:
            self.genome_stats[gen]["unique_reads"] += 1
        else:
            self.genome_stats[gen]["ambiguous_reads"] += 1
        return None

    def __dump_summary(self) -> None:
        """counts how many reads are unmapped, ambiguously mapped and uniquly mapped"""
        for read,mapping in self.summary.items():
            if mapping[0] is True:
                self.reads_stats["unique_mapped_reads"] += 1
                self.__incriment_genome(str(mapping[1][0]), True)
            else:
                if mapping[1] == []:
                    self.reads_stats["unmapped_reads"] += 1
                else:
                    self.reads_stats["ambiguous_mapped_reads"] += 1
                    for gen in mapping[1]:
                        self.__incriment_genome(str(gen), False)
        return None

    def get_summaries(self) -> Tuple[Dict[str, int], Dict[str, Dict[str, int]]]:
        """getter function for the genome summary and reas summary"""
        self.__dump_summary()
        return self.reads_stats, self.genome_stats
