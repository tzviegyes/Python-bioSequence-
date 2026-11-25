Biosequence Pseudo-Alignment Tool

A command-line bioinformatics tool in Python for building k-mer reference databases
    and pseudo-aligning NGS reads against reference genomes.

What It Does
- Builds k-mer reference databases (.kdb) from FASTA (.fa) genome files, with optional similarity filtering.
- Runs pseudo-alignment on FASTQ NGS reads against a reference file or directly from a genome file.
- Supports quality thresholds, unique/ambiguous thresholds, max-genome caps, and reverse-complement orientation selection.
- Prints or exports JSON summaries for reference databases and alignment results.

Tools & Languages
- Python
- Standard library modules (argparse, sys, json, typing)

How to Run
1. Make sure numpy is installed
2. Find full instructions in: python main.py --help 
3. A few commands for example:
    python main.py -t dumpalign -g example_input/reference_example_for_display.fa -k 31 --reads example_input/small_example_se_k75_n1000.fq
    python main.py -t dumpalign -g example_input/reference_example_for_display.fa -k 31 --reads example_input/small_example_se_k75_n1000.fq --min-read-quality 29 --min-kmer-quality 31.5
    python main.py -t dumpalign -g example_input/reference_example_for_display.fa -k 7 --reads example_input/small_example_se_k75_n1000.fq --max-genomes 2



