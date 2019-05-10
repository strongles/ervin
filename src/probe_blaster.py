from Bio.Blast.Applications import NcbiblastnCommandline
from ervin_utils import format_timestamp_for_filename
from probe_data import ProbeData
import progressbar
import argparse
import os

from exceptions import InvalidPathException

TEMP_FASTA_FILE = "/tmp/temp.fasta"
TEMP_TBLASTN_OUTPUT = "/tmp/temp_tblastn.tsv"

DEFAULT_OUTPUT_DIR = os.path.join(os.getcwd(), "OUTPUT")

DEFAULT_ALIGNMENT_LENGTH_THRESHOLD = 400
DEFAULT_E_VALUE_THRESHOLD = 0.009


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        help="Source fasta file containing the sample "
                             "probe records to run through tblastn",
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-db', "--database_path",
                        help="Path to the genome database against which "
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


def read_and_sanitise_raw_data(filename):
    with open(filename) as fasta_input:
        raw_data = fasta_input.readlines()
    stripped_newlines_data = [line.strip() for line in raw_data]
    purged_empty_lines_data = [line for line in stripped_newlines_data if line != ""]
    return purged_empty_lines_data


def read_probe_records_from_file(filename):
    input_data = read_and_sanitise_raw_data(filename)
    fasta_titles = [line for line in input_data if ">" in line]
    title_indices = [input_data.index(title) for title in fasta_titles]
    parsed_file_data = []
    for i in range(len(title_indices)):
        try:
            sequence = "".join(input_data[title_indices[i]+1:title_indices[i+1]])
            record = {"title": input_data[title_indices[i]], "seq": sequence}
            parsed_file_data.append(record)
        except IndexError:
            sequence = "".join(input_data[title_indices[i] + 1:])
            record = {"title": input_data[title_indices[i]], "seq": sequence}
            parsed_file_data.append(record)
    return parsed_file_data


def print_probe_to_temp_fasta_file(probe):
    with open(TEMP_FASTA_FILE, 'w') as fasta_outfile:
        fasta_outfile.write(f"{probe['title']}\n{probe['seq']}\n")


def run_blast(probe, db, e_value_threshold):
    print_probe_to_temp_fasta_file(probe)
    command = NcbiblastnCommandline(cmd="tblastn",
                                    out=TEMP_TBLASTN_OUTPUT,
                                    outfmt="\"6 qseqid sseqid slen sstart send evalue length qseq sseq sframe\"",
                                    query=TEMP_FASTA_FILE,
                                    db=db,
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


def print_results(result_list, title, run_time, output_dir):
    if output_dir is None:
        output_dir = DEFAULT_OUTPUT_DIR
    if not os.path.isdir(output_dir):
        if os.path.exists(output_dir):
            raise InvalidPathException(f"Invalid output path provided: {output_dir}")
        else:
            os.makedirs(output_dir)
    output_filename = f"{title}_{run_time}.tsv"
    output_filepath = os.path.join(output_dir, output_filename)
    with open(output_filepath, 'w') as result_out:
        for result in sorted(result_list):
            result_out.write(result.to_tsv())


def run():
    args = parse_args()
    probe_records = read_probe_records_from_file(args.file.name)
    run_time = format_timestamp_for_filename()
    with progressbar.ProgressBar(max_value=len(probe_records), type="percentage") as bar:
        for count, probe_record in enumerate(probe_records):
            blast_result = run_blast(probe_record, args.database_path, args.e_value)
            filtered_results = filter_results(blast_result, args)
            print_results(filtered_results, probe_record["title"], run_time, args.output_dir)
            bar.update(count)


if __name__ == "__main__":
    run()
