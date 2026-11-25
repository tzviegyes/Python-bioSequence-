import pytest
from Reference import Reference
from NGSReads import NgsReads
from Kmer_reference import KmerReference
from PsuedoAlignment import PsuedoAlginment
from main import Main
import os
import tempfile

class TestsProject():
    def test_invalid_input(self):
        assert(Reference(['organ_1'], "CTGTCTCTG", 3)._valid_all_input(), False) # name not string
        assert(Reference("organ_2", "CTGTCL", 3)._valid_all_input(), False) # non ACGTN letter
        assert(Reference("organ_3", "CTGTCTHT", "m")._valid_all_input(), False) # non integer K
        assert(Reference("organ_4", ["a", "b", "c", 'd'], 3)._valid_all_input(), False) # non string sequence
        assert(Reference("organ_5", "CTGT", 5)._valid_all_input(), False) # K > length
        assert(Reference("organ_6", "CTGTGG", 5)._valid_all_input(), True) # OK

    def test_ref_kmer_dict(self):
        genome_1 = Reference("genome_1", "ACTGGGGACT", 3)
        assert(genome_1.ref_kmer_dict,{'ACT': [0, 7], 'CTG': [1], 'TGG': [2], 'GGG': [3, 4], 'GGA': [5], 'GAC': [6]})
        genome_2 = Reference("genome_2", "NACTG", 3)
        assert(genome_2.ref_kmer_dict,{'ACT': [1], 'CTG': [2]})
        genome_3 = Reference("genome_3", "NACANTAG", 1)
        assert(genome_3.ref_kmer_dict,{'A': [1, 3, 6], 'C': [2], 'T': [5], 'G': [7]})
        genome_4 = Reference("genome_4", "AGNACANTAG", 3)
        assert(genome_4.ref_kmer_dict,{'ACA': [3], 'TAG': [7]})
        sequence = "ATAGGGATCGATAGGGATCGATAGGGATCGNNNATAGGGATCGATAGGGATCGATAGGGATCCNNNNNNNATAGGGATCGATAGGGATCGATAGGGATCG"
        genome_5 = Reference("genome", sequence, 30)
        assert(genome_5.ref_kmer_dict, {'ATAGGGATCGATAGGGATCGATAGGGATCG': [0, 70], 'ATAGGGATCGATAGGGATCGATAGGGATCC': [33]})

    def test_algorithm_1(self):
        genome1_1 = Reference('genome1', "GGGGGG", 4)
        genome1_2 = Reference('genome2', "AAATTT", 4)
        kmer_reference1 = KmerReference([genome1_1, genome1_2], 4)
        ngs_read1 = [NgsReads("name", "TAGGCAT")]
        algo1 = PsuedoAlginment(ngs_read1, kmer_reference1)
        mapped = algo1.get_mapping_for_read()
        assert((False, []), mapped)

    def test_algorithm_2(self):
        genome2_1 = Reference('genome1', "ATGGCTA", 4)
        genome2_2 = Reference('genome2', "CTAT", 4)
        kmer_reference2 = KmerReference([genome2_1, genome2_2], 4)
        ngs_read1 = [NgsReads("name","ATGGCTAT")]
        algo2 = PsuedoAlginment(ngs_read1, kmer_reference2)
        mapped = algo2.get_mapping_for_read()
        assert(('True', 'genome1'), (str(mapped[0]), str(mapped[1][0])))

    def test_algorithm_4(self):
        genome4_1 = Reference('genome1', "ATGCCNCGGGG", 4)
        genome4_2 = Reference('genome2', "GCCGNCGGGGNGCTAANGGGCT", 4)
        genome4_3 = Reference('genome3', "CCGG", 4)
        genome4_4 = Reference('genome4', "GGGCT", 4)
        genome4_5 = Reference('genome5', "GCTAA", 4)

        kmer_reference4 = KmerReference([genome4_1, genome4_2, genome4_3, genome4_4, genome4_5], 4)
        ngs_read4 = [NgsReads("name", "ATGCCGGGGCTAA")]
        algo4 = PsuedoAlginment(ngs_read4, kmer_reference4)
        mapped = algo4.get_mapping_for_read()
        ambiguouse_genomes = algo4.get_genome_stats()
        assert(False, mapped[0])
        assert(('genome1', 'genome2') ,(str(mapped[1][0]), str(mapped[1][1])))
        assert(2, len(mapped[1]))
        assert((1, 1), (ambiguouse_genomes[str(genome4_1)]["ambiguous_reads"], ambiguouse_genomes[str(genome4_2)]["ambiguous_reads"]))
        assert((0, 0, 0), (ambiguouse_genomes[str(genome4_3)]["ambiguous_reads"], ambiguouse_genomes[str(genome4_4)]["ambiguous_reads"], ambiguouse_genomes[str(genome4_5)]["ambiguous_reads"]))

    def test_skip_N_referenc(self):
        genome4_2 = Reference('genome2', "GCCGNCGGGGNGCTAANGGGCT", 4)
        assert(['GCCG', 'CGGG', 'GGGG', 'GCTA', 'CTAA', 'GGGC', 'GGCT'],
                         list(genome4_2.get_reference_kmer_dict().keys()))

    def test_read_quality_check(self):
        read1 = NgsReads('read1', 'ATGCAT', "]]]]]]")
        read2 = NgsReads('read2', 'ATGCAT', "5?IS]g")
        assert(60.0 ,read1.get_read_qulity())
        assert(45.0, read2.get_read_qulity())
        assert(['TGCA', 'GCAT'], read2.get_all_kmers(4, 40))

    def test_counters_for_qulity_check(self):
        gen1 = Reference('gen1', 'ATGCAT', 4)
        gen2 = Reference('gen2', 'GGGCAT', 4)
        gen3 = Reference('gen3', 'TAGCAT', 4)

        read1 = NgsReads('read1', 'ATGCAT', "++++++")
        read2 = NgsReads('read2', 'ATGCAT', "5?IS]g")
        kmer_reference1 = KmerReference([gen1, gen2, gen3], 4)
        ngs_read1 = [read1, read2]
        algo1 = PsuedoAlginment(ngs_read1, kmer_reference1, 1, 1, 30, 40, 2)
        algo1.psuedo_algin_reads_list()

        assert(1, algo1.get_skiped_reads())
        assert(1, algo1.get_skipped_low_quality_kmer())
        assert(1, algo1.get_skipped_kmers_high_mapping())

    def test_reverse_complement1(self):
        read1 = NgsReads('read1', "ATCGCCAT")
        genome1 = Reference('gen1', 'ATCG', 4)
        genome2 = Reference('gen2', 'ATGGCG', 4)
        kmer_ref1 = KmerReference([genome1, genome2], 4)
        algo1 = PsuedoAlginment([read1], kmer_ref1, check_reverse=True).get_mapping_for_read()
        assert(('True', 'gen2'), (str(algo1[0]), str(algo1[1][0])))

    def test_tasks12347(self):
        # write temporary genome file
        with (tempfile.NamedTemporaryFile(mode='w+', delete=False) as genomes1,
              tempfile.NamedTemporaryFile('w+', delete=False) as reads4):
            genomes_1_name = genomes1.name
            genomes1.write(">genome1\nGGGGGG\n")
            genomes1.write(">genome2\nAAATTT\n")
            genomes1.seek(0)

            reads_4_name = reads4.name
            reads4.write("@name\nAAATT\n")
            reads4.write("+\n@#$%^\n")
            reads4.seek(0)

        try:
            # Create another temporary file for the output
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_output:
                refrence_1_name = temp_output.name

            runner = Main({'k': 4})
            # test function 1. writes output to temp_output
            runner._building_a_reference_1(genomes_1_name, 4, refrence_1_name)

            # test function 2 from output of task 1
            output_string_task_2 = runner._summeries_reference_2(refrence_1_name)
            printable_output_2 = output_string_task_2.replace('\n', '').split('},')
            assert('{"Kmers":{"GGGG":{"genome1":[0,1,2]', printable_output_2[0])
            assert('"Summary":{"genome1":{"total_bases":6,"unique_Kmers":1,"multi_mapping_kmers":0', printable_output_2[4])

            # test 3 from same genome file
            printable_output_3 = runner._build_dump_reference_3(genomes_1_name, 4).replace('\n', '').split('},')
            assert('{"Kmers":{"GGGG":{"genome1":[0,1,2]', printable_output_3[0])
            assert('"Summary":{"genome1":{"total_bases":6,"unique_Kmers":1,"multi_mapping_kmers":0', printable_output_3[4])
            # test 4 from output of 1
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_output:
                align_4_name = temp_output.name
                runner._psuedo_align_pre_built_4(refrence_1_name, align_4_name, reads_4_name, 1, 1, 256, 0, -1, False)
                # test 6 from output of 4
                result = runner._dump_pre_built_align_6(align_4_name).replace('\n', '')
                assert('{"Statistics":{"unique_mapped_reads":1,"ambiguous_mapped_reads":0,"unmapped_reads":0}', result[:85])
                assert('"Summary":{"genome2":{"unique_reads":1,"ambiguous_reads":0}}}', result[86:])

            # test 7 with refrence file from function 1
            with open(refrence_1_name, "r"):
                result7 = runner._algorithm_prebuilt_dump_7(refrence_1_name, reads_4_name, 1, 1, 256, 0, -1,
                                                             False).replace('\n', '')
                assert(
                    '{"Statistics":{"unique_mapped_reads":1,"ambiguous_mapped_reads":0,"unmapped_reads":0}', result7[:85])
                assert('"Summary":{"genome1":{"unique_reads":0,"ambiguous_reads":0},'
                                 '"genome2":{"unique_reads":1,"ambiguous_reads":0}}}', result7[86:])
        finally:
            # Clean up temporary files
            os.remove(genomes_1_name)
            os.remove(refrence_1_name)
            os.remove(reads_4_name)
            os.remove(align_4_name)

    def test_tasks56(self):
        # write temporary genome file
        with (tempfile.NamedTemporaryFile('w+', delete=False) as genomes5,
              tempfile.NamedTemporaryFile('w+', delete=False) as reads5):
            genomes_5_name = genomes5.name
            genomes5.write(">genome1\nGGGGGG\n")
            genomes5.write(">genome2\nAAATTT\n")
            genomes5.seek(0)

            reads_5_name = reads5.name
            reads5.write("@name\nAAATT\n")
            reads5.write("+\n@#$%^&\n")
            reads5.seek(0)
        try:
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_output:
                align5 = temp_output.name
                runner = Main({'k': 4})
                runner._psuedo_align_to_file_5(genomes_5_name, 4, align5, reads_5_name, 1, 1, 256, 0, -1, False)
                # test 6 from output of 5
                result = runner._dump_pre_built_align_6(align5).replace('\n', '')
                assert(
                    '{"Statistics":{"unique_mapped_reads":1,"ambiguous_mapped_reads":0,"unmapped_reads":0}',
                    result[:85])
                assert('"Summary":{"genome2":{"unique_reads":1,"ambiguous_reads":0}}}', result[86:])
        finally:
            os.remove(genomes_5_name)
            os.remove(reads_5_name)
            os.remove(align5)

    def test_task8(self):
        # write temporary genome file
        with (tempfile.NamedTemporaryFile('w+', delete=False) as genomes8,
              tempfile.NamedTemporaryFile('w+', delete=False) as reads8):
            genomes_8_name = genomes8.name
            genomes8.write(">genome1\nGGGGGG\n")
            genomes8.write(">genome2\nAAATTT\n")
            genomes8.seek(0)

            reads_8_name = reads8.name
            reads8.write("@name\nAAATT\n")
            reads8.write("+\n@#$%^&\n")
            reads8.seek(0)
        try:
            runner = Main({'k': 4})
            result8 = runner._build_align_dump_8(genomes_8_name, 4, reads_8_name, 1, 1, 256, 0, -1,
                                                  False).replace('\n', '')
            assert(
                '{"Statistics":{"unique_mapped_reads":1,"ambiguous_mapped_reads":0,"unmapped_reads":0}', result8[:85])
            assert('"Summary":{"genome1":{"unique_reads":0,"ambiguous_reads":0},'
                             '"genome2":{"unique_reads":1,"ambiguous_reads":0}}}', result8[86:])
        finally:
            os.remove(genomes_8_name)
            os.remove(reads_8_name)
            # write temporary genome file
            with (tempfile.NamedTemporaryFile('w+', delete=False) as genomes5,
                  tempfile.NamedTemporaryFile('w+', delete=False) as reads5):
                genomes_5_name = genomes5.name
                genomes5.write(">genome1\nGGGGGG\n")
                genomes5.write(">genome2\nAAATTT\n")
                genomes5.seek(0)

                reads_5_name = reads5.name
                reads5.write("@name\nAAATT\n")
                reads5.write("+\n@#$%^&\n")
                reads5.seek(0)
            try:
                runner = Main({'k': 4})
                result = runner._build_align_dump_8(genomes_5_name, 4, reads_5_name, 1, 1, 256, 0, -1, False)
                result = result.replace('\n', '')
                assert(
                    '{"Statistics":{"unique_mapped_reads":1,"ambiguous_mapped_reads":0,"unmapped_reads":0}',
                    result[:85])
                assert('Summary":{"genome1":{"unique_reads":0,"ambiguous_reads":0},"genome2":'
                                 '{"unique_reads":1,"ambiguous_reads":0}}}', result[87:])
            finally:
                os.remove(genomes_5_name)
                os.remove(reads_5_name)

    def test_filter_similar(self):
        genome1 = Reference("gen1", 'ATGCATG', 4)
        genome2 = Reference("gen2", 'ATGCA', 4)
        genome3 = Reference('gen3', 'TATGC', 4)

        kmer_ref1 = KmerReference([genome1, genome2, genome3], 4, 0.45)
        gen_filter_stats = kmer_ref1.get_filter_stats()
        assert(gen_filter_stats['gen1']['kept'], 'yes')
        assert(gen_filter_stats['gen2']['kept'], 'no')
        assert(gen_filter_stats['gen3']['kept'], 'no')
        assert(1, len(kmer_ref1.get_genomes()))

        kmer_ref2 = KmerReference([genome1, genome2, genome3], 4, 0.95)
        gen_filter_stats = kmer_ref2.get_filter_stats()
        assert(gen_filter_stats['gen1']['kept'], 'yes')
        assert(gen_filter_stats['gen2']['kept'], 'no')
        assert(gen_filter_stats['gen3']['kept'], 'yes')
        assert(2, len(kmer_ref2.get_genomes()))

if __name__ == '__main__':
    TestsProject()