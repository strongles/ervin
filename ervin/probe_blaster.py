from .ervin_utils import DEFAULT_OUTPUT_DIR
from .ervin_utils import TEMP_FASTA_FILE
from .ervin_utils import TEMP_PROBE_BLASTER
from .ervin_utils import TEMP_TBLASTN_OUTPUT
from .ervin_utils import get_config
from .ervin_utils import format_timestamp_for_filename
from .ervin_utils import read_from_fasta_file
from .ervin_utils import total_result_records

from .exceptions import InvalidPathException

from .probe_data import ProbeData

from Bio.Blast.Applications import NcbiblastnCommandline
from collections import namedtuple

import argparse
import os
import progressbar


DEFAULT_ALIGNMENT_LENGTH_THRESHOLD = 400
DEFAULT_E_VALUE_THRESHOLD = 0.009

Args = namedtuple("Args", "filename output_dir alignment_len_threshold e_value")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        help="Source fasta file containing the sample "
                             "probe records to run through tblastn",
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-gdb', "--genome_database",
                        help="Name of the genome database against which "
                             "the probe records are to be BLASTed",
                        type=str,
                        required=True)
    parser.add_argument("-o", "--output_dir",
                        help="Location to which to write the result files",
                        type=str,
                        required=False)
    parser.add_argument("-a", "--alignment_len_threshold",
                        help="Minimum length threshold that BLAST result "
                             "alignment sequence lengths should exceed",
                        type=int,
                        required=False,
                        default=DEFAULT_ALIGNMENT_LENGTH_THRESHOLD)
    parser.add_argument("-e", "--e_value",
                        help="Maximum e-value by which to threshold the "
                             "results returned by the BLAST run",
                        type=float,
                        required=False,
                        default=DEFAULT_E_VALUE_THRESHOLD)
    return parser.parse_args()


def print_probe_to_temp_fasta_file(probe):
    with open(TEMP_FASTA_FILE, 'w') as fasta_outfile:
        fasta_outfile.write(f"{probe['title']}\n{probe['seq']}\n")


def run_blast(probe, db, e_value_threshold):
    config = get_config()
    print_probe_to_temp_fasta_file(probe)
    command = NcbiblastnCommandline(cmd="tblastn",
                                    out=TEMP_TBLASTN_OUTPUT,
                                    outfmt="\"6 qseqid sseqid slen sstart send "
                                           "evalue length qseq sseq sframe\"",
                                    query=TEMP_FASTA_FILE,
                                    db=f"{config.genome_db_storage}{db}",
                                    evalue=e_value_threshold)
    command()
    probe_records = []
    with open(TEMP_TBLASTN_OUTPUT) as blast_output:
        for line in blast_output:
            probe_records.append(ProbeData(line.strip()))
    return probe_records


def _length_requirement(probe_list, args):
    return [probe for probe in probe_list if probe.alignment_length > args.alignment_len_threshold]


def filter_results(hits, args):
    filters = [
        _length_requirement,
    ]
    for filter_case in filters:
        hits = filter_case(hits, args)
    return hits


def print_results(result_list, title, run_time, output_dir=DEFAULT_OUTPUT_DIR):
    if not os.path.isdir(output_dir):
        if os.path.exists(output_dir):
            raise InvalidPathException(f"Invalid output path provided: {output_dir}")
        else:
            os.makedirs(output_dir)
    output_filename = f"{title.replace('>', '')}_{run_time}.tsv"
    output_filepath = os.path.join(output_dir, output_filename)
    with open(output_filepath, 'w') as result_out:
        for result in sorted(result_list):
            result_out.write(result.to_tsv())

    return output_filepath


def run_probe_blaster(file, genome_db, align_threshold, e_value, run_ts=None, output_dir=TEMP_PROBE_BLASTER):
    run_time = run_ts if run_ts else format_timestamp_for_filename()
    align_len = align_threshold if align_threshold else DEFAULT_ALIGNMENT_LENGTH_THRESHOLD
    e_val_threshold = e_value if e_value else DEFAULT_E_VALUE_THRESHOLD
    args = Args(file, output_dir, align_len, e_val_threshold)
    probe_records = read_from_fasta_file(file)
    output_filepaths = []
    with progressbar.ProgressBar(max_value=len(probe_records),
                                 type="percentage",
                                 prefix="Blasting against Genome DB: ") as bar:
        for count, probe_record in enumerate(probe_records):
            bar.update(count)
            blast_result = run_blast(probe_record, genome_db, e_val_threshold)
            filtered_results = filter_results(blast_result, args)
            output_file = print_results(filtered_results,
                                        probe_record["title"],
                                        run_time, output_dir)
            output_filepaths.append(output_file)
    return output_filepaths, total_result_records(output_filepaths)


if __name__ == "__main__":
    args = parse_args()
    run_probe_blaster(args.file.name, args.genome_database, args.alignment_len_threshold,
                      args.e_value, args.output_dir)
