import argparse
import sys
import FileHandling
from Kmer_reference import KmerReference
from PsuedoAlignment import PsuedoAlginment, AlignmentData
from typing import Optional
import json

PROGRAM_DESCRIPTION = ("This program has the following tools:\n"
                 "-t reference - Build a K-mer reference database,\n"
                 "-t dumpref - Dump a K-mer reference database.\n"
                 "-t align - Preform psuedo alginment for a NGS reads file.\n"
                 "-t dumpalign - Dump the output of the psuedo alginment.\n\n"
                 "reference: Builds a K-mer reference database.\n"
                 "    Enter -g <genome_file_name.fa> -k <K-mer size> -r <output_file.kdb>\n\n"
                 "dumpref: Produces a flat summary of K-mer reference file.\n"
                 "    For dumping a refrence file enter -r <reference_file.kdb>\n"
                 "    For creating a K-mer reference database and dumping it - enter  -g <genome_file_name.fa> -k <K-mer size>\n\n"
                 "  For filtering similar genomes enter --filter-similar and --similarity-threshold THRESHOLD.\n\n"
                 "align: Runs the Psuedo Alignment algorithm.\n"
                 "    For building a reference database and aligning - enter -g <genome_file_name.fa> -k <K-mer size>\n"
                 "      -a <alignment_output.aln> --reads <NGS_reads.fastq> optional: [-m UNIQE_THRESHOLD] [-p AMBIGUOUS_THRESHOLD]\n"
                 "    For aligning a provided reference file - enter -r <reference_file.kdb> -a <alignment_output.aln>\n"
                 "      --reads <NGS_reads.fastq> optional: [-m UNIQE_THRESHOLD] [-p AMBIGUOUS_THRESHOLD].\n\n"
                 "dumpalign: Dumps a summery of Psuedo Alignment algorithm.\n"
                 "    To produce a flat summary of a given alignment file - enter -a <alignment_file.aln>\n"
                 "    To align a given reference file and summereize result - enter -r <reference_file.kdb>\n"
                 "      --reads <NGS_reads.fastq> optional: [-m UNIQE_THRESHOLD] [-p AMBIGUOUS_THRESHOLD].\n"
                 "    To build a reference database, align and dump result - enter -g <genome_file_name.fa> -k <K-mer size>\n"
                 "      --reads <NGS_reads.fastq> optional: [-m UNIQE_THRESHOLD] [-p AMBIGUOUS_THRESHOLD].\n\n"
                 "  For align and dumpalign - enter --min-read-quality, --min-kmer-quality, --max-genomes for reads quality filtering.\n"
                 "    eter --reverse-complement for alignment with the direction of read that gives better result.\n"
                 "    enter --help for more information.\n")
def parser() -> argparse.Namespace:
    """Parses the arguments from user's input"""
    parser = argparse.ArgumentParser(description= PROGRAM_DESCRIPTION,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     add_help=False, prog='Biosequence project',)
    parser.add_argument('-h', '--help', action="help", help="Use the following tools to analyse DNA sequences:")
    parser.add_argument('-g', '--genomefile', type=str, help="A .fa file of reference genomes.")
    parser.add_argument('-k', '--kmer-size', type=int, help="Lenght of K-mers to search.")
    parser.add_argument('-r', '--referencefile', type=str, help="A .kdb file of the K-mer reference database.")
    parser.add_argument('-t', '--task', choices=['dumpref', 'align', 'dumpalign', 'reference'],
                        required=True,help="Enter --help for more information about the tools.\n")
    parser.add_argument('-a', '--alignfile', type=str, help="A file to dump the Psuedo Alignment's summery.")
    parser.add_argument('--reads', type=str, help="A FASTQ file of NGS reads.")
    parser.add_argument('-m', '--unique-threshold', default=1, type=int, help="The unique threshold. default is 1")
    parser.add_argument('-p','--ambiguous-threhold', default=1, type=int, help="The amgiguouse threshold. default is 1")
    parser.add_argument('--min-read-quality', default=256, type=float, help="The minimum threshold for the quality of the NGS DNA read")
    parser.add_argument('--min-kmer-quality', default=0, type=float, help='The minimum threshold for the quality of a K-mer')
    parser.add_argument('--max-genomes', default=-1, type=int, help='The maximum number of genomes a K-mer can be mapped to,\n'
                                                                    'otherwise the algorithm wiil ignor the K-mer\n')
    parser.add_argument('--reverse-complement', action='store_true',
                        help="enter flag if you want program to choose the orientation with higher specific Kmer count\n")
    parser.add_argument('--filter-similar', action='store_true',
                        help="enter flag if you want program to filter similar genome refrences.\n"
                             "Enter --similarity-threshold THRESHOLD for the simularity threshold score to filter.\n"
                             "default is 0.95\n")
    parser.add_argument('--similarity-threshold', default=0.95, type=float,
                        help= "Choose threshold for filtering. Enter a float larger than 0 and smaller than 1\n"
                              "default is 0.95\n")
    return parser.parse_args()

class Main:
    """The class that runs the program"""
    def __init__(self, arguments):
        """constructor for Main object, recives the input arguments"""
        self.arguments = arguments

    def _building_a_reference_1(self, genome_file: str, kmer: int, reference_output: str, if_to_filter=False, filter_threshold=0.95) -> None:
        """Preformes task 1. Builds a reference file frome genome file and K-mer length"""
        # print("Extracting genome file")
        genome_list = FileHandling.make_list_of_genomes(genome_file, kmer)
        # print("Building a reference")
        if not if_to_filter:
            kmer_reference = KmerReference(genome_list, kmer)
        else:
            kmer_reference = KmerReference(genome_list, kmer, filter_threshold)
            # print(json.dumps({"Similarity": kmer_reference.get_filter_stats()}, indent=4, separators=(",", ":")))
        filter_stats = kmer_reference.get_filter_stats()
        # print("Writing reference file")
        FileHandling.write_data_to_file(reference_output,
                {"reference_database":kmer_reference.get_database_for_printing(), "filter_stats": filter_stats})
        # print("Task completed successfully\n")
        return None

    def _summeries_reference_2(self, refrence_file: str) -> Optional[str]:
        """Preforms task 2. dumps a reference file"""
        # print("Extracting data from file")
        data = FileHandling.read_data_file(refrence_file)
        # print("Processing the data")
        kmer_reference = KmerReference(data)
        kmer_reference_for_print = kmer_reference.get_database_for_printing()
        summary = {"Kmers": kmer_reference_for_print, "Summary": kmer_reference.get_genome_summery()}
        filter_stats = kmer_reference.get_filter_stats()
        if filter_stats != {}:
            summary.update({"Similarity": filter_stats})
        if not __name__ == "__main__": # for testing
            return json.dumps(summary, indent=0, separators=(",", ":"))
        else:
            print(json.dumps(summary, indent=4)) # , separators=(",", ":")
            # print("Task completed successfully\n")
        return None

    def _build_dump_reference_3(self, genome_file: str, kmer: int, if_to_filter=False, filter_threshold=0.95) -> Optional[str]:
        """Preforms task 3. Builds a reference from genome fa file and K-mer lenght.
        Dumps the output."""
        # print("Extracting genome file")
        genome_list = FileHandling.make_list_of_genomes(genome_file, kmer)
        # print("Building a reference")
        if not if_to_filter:
            kmer_reference = KmerReference(genome_list, kmer)
        else:
            kmer_reference = KmerReference(genome_list, kmer, filter_threshold)
        kmer_reference_for_print = kmer_reference.get_database_for_printing()
        summary = {"Kmers": kmer_reference_for_print, "Summary": kmer_reference.get_genome_summery()}
        if not __name__ == "__main__": # for testing
            return json.dumps(summary, indent=0, separators=(",", ":"))
        else:
            if if_to_filter:
                summary.update({"Similarity": kmer_reference.get_filter_stats()})
            print(json.dumps(summary, indent=4))
        return None

    def _psuedo_align_pre_built_4(self, reference_file: str, align_file: str, reads_file: str, *args) -> None:
        """Preforms task 4. Psuedo aligns a refernce file with a fq NGS reads file"""
        data = FileHandling.read_data_file(reference_file)
        kmer_reference = KmerReference(data)
        reads = FileHandling.make_list_of_reads(reads_file)
        mapping = PsuedoAlginment(reads, kmer_reference, *args)
        FileHandling.write_data_to_file(align_file, mapping.get_psuedo_algin_reads_list())
        return None

    def _psuedo_align_to_file_5(self, genome_file: str, kmer: int, align_file: str, reads_file: str, *args) -> None:
        """Prefoems task 5. Writes a psuedo alignment file from given fa genome file, K-mer length  and fq NGS reads file"""
        genome_list = FileHandling.make_list_of_genomes(genome_file, kmer)
        kmer_reference = KmerReference(genome_list, kmer)
        reads = FileHandling.make_list_of_reads(reads_file)
        mapping = PsuedoAlginment(reads, kmer_reference, *args)
        FileHandling.write_data_to_file(align_file, mapping.get_psuedo_algin_reads_list())
        return None

    def _dump_pre_built_align_6(self, align_file: str) -> Optional[str]:
        """Preforms task 6. Dumps a psuedo alignment summary file"""
        # print("Extracting alignment file")
        data = FileHandling.read_data_file(align_file)
        # print("Proccesing data")
        summary = AlignmentData(data)
        reads_summary, genome_summary = summary.get_summaries()
        summary_reads_genomes = {"Statistics": reads_summary, "Summary": genome_summary}
        if not __name__ == "__main__":
            return json.dumps(summary_reads_genomes, indent=0, separators=(",", ":"))
        else:
            print(json.dumps(summary_reads_genomes, indent=4))
            # print("Task completed successfully\n")
            return None

    def _algorithm_prebuilt_dump_7(self, reference_file: str, reads_file: str, *args) -> Optional[str]:
        """Preforms task 7. Recives a reference file, and fq NGS reads file.
        Psuedo aligns the genomes and reads then prints the alignmrnt summary"""
        data = FileHandling.read_data_file(reference_file)
        kmer_reference = KmerReference(data)
        load_quality = True if args[2] != 256 else False
        reads = FileHandling.make_list_of_reads(reads_file, load_quality)
        mapping = PsuedoAlginment(reads, kmer_reference, *args)
        summary = {"Statistics": mapping.get_reads_stats(), "Summary": mapping.get_genome_stats()}
        filter_stats = kmer_reference.get_filter_stats()
        if filter_stats != {}:
            summary.update({"Similarity": filter_stats})
        if not __name__ == '__main__':
            return json.dumps(summary, indent=0, separators=(",", ":"))
        else:
            print(json.dumps(summary, indent=4))
        return None

    def _build_align_dump_8(self, genome_file: str, kmer: int, reads_file: str, *args) -> Optional[str]:
        """Preforms task 8. Builds a reference frome genome file and K-mer lenght.
        Psuedo aligns with reads file and prints the alignmrnt summary"""
        genome_list = FileHandling.make_list_of_genomes(genome_file, kmer)
        kmer_reference = KmerReference(genome_list, kmer)
        load_quality = True if args[2] != 256 else False
        reads = FileHandling.make_list_of_reads(reads_file, load_quality)
        mapping = PsuedoAlginment(reads, kmer_reference, *args)
        summary = {"Statistics": mapping.get_reads_stats(), "Summary": mapping.get_genome_stats()}
        if not __name__ == "__main__":
            return json.dumps(summary, indent=0, separators=(",", ":"))
        else:
            print(json.dumps(summary, indent=4))
        return None

    def __handle_t_reference(self) -> None: # task 1
        """Calls the right tool function if user input is -t reference"""
        genome_file = self.arguments.genomefile
        kmer_size = self.arguments.kmer_size
        ref_file = self.arguments.referencefile
        similar_threshold = self.arguments.similarity_threshold
        if ref_file is None:
            print("\nFor dumping existiing reference - enter the filename or path.\n"
                  "For creating new reference - enter output file name, genome file, and K-mer size\n")
            exit(1)
        else:
            if genome_file is None or kmer_size is None:
                print(
                    "\nWhen using reference tool, input K-mer size, reference file or path, and output file or path.\n"
                    "see --help for more information.\n")
                exit(1)
            else:
                if not 0 < similar_threshold < 1 or not isinstance(similar_threshold, float):
                    print("Threshold must be a float larger than 0 and smaller than 1.\n")
                    exit(1)
                else:
                    self._building_a_reference_1(genome_file, kmer_size, ref_file, self.arguments.filter_similar,
                                                 similar_threshold)

    def __handle_t_dumpref(self) -> None: # tasks 2,3
        """Calls the right tool function if user input is -t dumpref"""
        similar_threshold = self.arguments.similarity_threshold
        if self.arguments.referencefile is not None:
            self._summeries_reference_2(self.arguments.referencefile)

        elif self.arguments.genomefile is None or self.arguments.kmer_size is None:
            print("\nFor dumping a referene - enter a genome .fa file and a K-mer lenght\n")
            sys.exit(1)
        else:
            if not 0 < similar_threshold < 1 or not isinstance(similar_threshold, float):
                print("Threshold must be a float larger than 0 and smaller than 1.\n")
                exit(1)
            else:
                self._build_dump_reference_3(self.arguments.genomefile, self.arguments.kmer_size, self.arguments.filter_similar,
                                             similar_threshold)

    def __handle_t_align(self) -> None: # tasks 4,5
        """Calls the right tool function if user input is -t align"""
        psuedo_align_argumeents = [self.arguments.unique_threshold, self.arguments.ambiguous_threhold, self.arguments.min_read_quality,
                                   self.arguments.min_kmer_quality, self.arguments.max_genomes, self.arguments.reverse_complement]
        if self.arguments.reads is None or self.arguments.alignfile is None:
            print("\nFor align tool - enter align file name or path and NGS reads file or path\n")
            sys.exit(1)
        if self.arguments.referencefile is not None:
            self._psuedo_align_pre_built_4(self.arguments.referencefile, self.arguments.alignfile, self.arguments.reads,
                                           *psuedo_align_argumeents)
        else:
            if self.arguments.genomefile is None or self.arguments.kmer_size is None:
                print(
                    "\nFor align tool - enter genome file or path and K-mer lenght, or enter reference file or path\n")
                sys.exit(1)
            else:
                self._psuedo_align_to_file_5(self.arguments.genomefile, self.arguments.kmer_size, self.arguments.alignfile, self.arguments.reads,
                                             *psuedo_align_argumeents)

    def __handle_t_dumpalign(self) -> None: # tasks 6,7,8
        """Calls the right tool function if user input is -t dumpalign"""
        if self.arguments.alignfile is not None:
            self._dump_pre_built_align_6(self.arguments.alignfile)
        elif self.arguments.reads is not None:
            psuedo_align_argumeents = [self.arguments.unique_threshold, self.arguments.ambiguous_threhold, self.arguments.min_read_quality,
                                       self.arguments.min_kmer_quality, self.arguments.max_genomes, self.arguments.reverse_complement]
            if self.arguments.referencefile is not None:
                self._algorithm_prebuilt_dump_7(self.arguments.referencefile, self.arguments.reads, *psuedo_align_argumeents)
            elif self.arguments.genomefile is not None and self.arguments.kmer_size is not None:
                self._build_align_dump_8(self.arguments.genomefile, self.arguments.kmer_size, self.arguments.reads,
                                         *psuedo_align_argumeents)
            else:
                print("\nFor dump alignment - enter alignment fiile or path,\n"
                      "or enter reference file and NGS reads file,\n"
                      "or enter genome file, K-mer length and NGS reads file.\n")
                sys.exit(1)
        else:
            print("\nFor dump alignment - enter alignment fiile or path,\n"
                  "or enter reference file and NGS reads file,\n"
                  "or enter genome file, K-mer length and NGS reads file.\n")
            sys.exit(1)

    def __asked_one_tool(self) -> bool:
        """Makes sur user chose exactly one tool"""
        if sys.argv.count("-t") > 1:
            print("\nYou can use only one tool in each command\n")
            return False
        if sys.argv.count("-t") == 0:
            print("\nYou didn't choose a tool. see --help for more information\n")
            return False
        return True

    def __valid_extra_quality(self) -> bool:
        """check if the extra quality check's input have valid value"""
        if self.arguments.min_read_quality < 0 or self.arguments.min_kmer_quality < 0:
            print("Quality can't be less than 0\n")
            return False
        if self.arguments.max_genomes < 1 and self.arguments.max_genomes != -1:
            print("Max genomes allowed needs to be at least 1\n")
            return False
        return True

    def run_program(self) -> None:
        """Identifies the tool the user chose and call the relavent function"""
        if not self.__asked_one_tool() or not self.__valid_extra_quality():
            exit(1)
        if self.arguments.task == 'reference':
            self.__handle_t_reference()
        if self.arguments.task == 'dumpref':
            self.__handle_t_dumpref()
        if self.arguments.task == 'align':
            self.__handle_t_align()
        if self.arguments.task == 'dumpalign':
            self.__handle_t_dumpalign()

if __name__ == '__main__':
    arguments = parser()
    Main(arguments).run_program()
