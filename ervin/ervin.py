from .probe_finder import run_probe_finder

from .probe_blaster import run_probe_blaster

from .ervin_utils import DEFAULT_OUTPUT_DIR
from .ervin_utils import ensure_output_dir_exists
from .ervin_utils import format_timestamp_for_filename

from .virus_blaster import VIRUS_DB_DEFAULT
from .virus_blaster import run_virus_blaster

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


def write_summary_file(probe_count, fasta_count, virus_to_counts, output_dir, run_ts):
    destination_dir = ensure_output_dir_exists(output_dir)
    summary_path = destination_dir / f"summary_{run_ts}.txt"
    with open(summary_path, "w") as summary:
        summary.write(f"PROBE HITS FOUND: {probe_count}\n")
        summary.write(f"PROBE HITS POST FILTERING: {fasta_count}\n")
        summary.write("VIRUS HITS BY NAME:\n")
        for virus_name, count in virus_to_counts.items():
            tab_char = "\t"
            summary.write(f"{tab_char}{virus_name}: {count}\n")


def run():
    args = get_args()
    run_ts = format_timestamp_for_filename()
    blasted_probes, probe_count = run_probe_blaster(args.file.name,
                                                    args.genome_database,
                                                    args.alignment_len_threshold,
                                                    args.e_value,
                                                    run_ts)
    probe_finder_fasta, fasta_count = run_probe_finder(blasted_probes,
                                                       args.alignment_len_threshold,
                                                       run_ts)
    _, virus_to_counts = run_virus_blaster(probe_finder_fasta,
                                           args.virus_database,
                                           args.output_dir,
                                           run_ts)
    write_summary_file(probe_count, fasta_count, virus_to_counts, args.output_dir, run_ts)
