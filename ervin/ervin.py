from .probe_finder import run_probe_finder

from .probe_blaster import run_probe_blaster

from .ervin_utils import DEFAULT_OUTPUT_DIR

from .virus_blaster import run_virus_blaster
from .virus_blaster import VIRUS_DB_DEFAULT

import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        help="Input probe file",
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument("-o", "--output_dir",
                        help="Location in which to write the output files",
                        type=str,
                        required=False, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("-gdb", "--genome_database",
                        help="Reference genome database",
                        type=str,
                        required=True)
    parser.add_argument("-vdb", "--virus_database",
                        help="Locally held Viruse reference database",
                        type=str,
                        required=False,
                        default=VIRUS_DB_DEFAULT)
    parser.add_argument("-a", "--alignment_len_threshold",
                        help="Minimum length threshold that BLAST result "
                             "alignment sequence lengths should exceed",
                        type=int,
                        required=False)
    parser.add_argument("-e", "--e_value",
                        help="Maximum e-value by which to threshold the "
                             "results returned by the BLAST run",
                        type=float,
                        required=False)
    return parser.parse_args()


def run():
    args = get_args()
    blasted_probes = run_probe_blaster(args.file.name,
                                       args.genome_database,
                                       args.alignment_len_threshold,
                                       args.e_value)
    probe_finder_fasta = run_probe_finder(blasted_probes, args.alignment_len_threshold)
    run_virus_blaster(probe_finder_fasta, args.virus_database, args.output_dir)
