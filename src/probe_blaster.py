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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        help="Source fasta file containing the sample probe records to run through tblastn",
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-db', "--database_path",
                        help="Path to the genome database against which the probe records are to be BLASTed",
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument("-o", "--output_dir",
                        help="Location to which to write the result files",
                        type=str,
                        required=False)
    return parser.parse_args()


def read_probe_records_from_file(filename):
    with open(filename) as fasta_input:
        raw_file_data = fasta_input.readlines()
    parsed_file_data = []
    for i in range(0, len(raw_file_data), 2):
        parsed_file_data.append({"title": raw_file_data[i].strip(), "seq": raw_file_data[i+1].strip()})
    return parsed_file_data


def print_probe_to_temp_fasta_file(probe):
    with open(TEMP_FASTA_FILE, 'w') as fasta_outfile:
        fasta_outfile.write(f"{probe['title']}\n{probe['seq']}\n")


def run_blast(probe, db):
    print_probe_to_temp_fasta_file(probe)
    command = NcbiblastnCommandline(cmd="tblastn",
                                    out=TEMP_TBLASTN_OUTPUT,
                                    outfmt="\"6 sseqid qseqid slen sstart send evalue length qseq sseq sframe\"",
                                    # outfmt=15,
                                    # outfmt=6,
                                    query=TEMP_FASTA_FILE,
                                    db=db,
                                    evalue=0.009)
    command()
    probe_records = []
    with open(TEMP_TBLASTN_OUTPUT) as blast_output:
        for line in blast_output:
            probe_records.append(ProbeData(line.strip()))
    return probe_records


def _length_requirement(probe_list):
    # filtered_hsps = [entry for entry in probe["hsps"] if entry["align_len"] > 400]
    # probe["hsps"] = filtered_hsps

    return [probe for probe in probe_list if probe.alignment_length > 400]


def filter_results(hits):
    filters = [
        _length_requirement,
    ]
    for filter_case in filters:
        hits = filter_case(hits)
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
        for result in result_list:
            result_out.write(result.to_tsv())


def run():
    args = parse_args()
    probe_records = read_probe_records_from_file(args.file.name)
    run_time = format_timestamp_for_filename()
    with progressbar.ProgressBar(max_value=len(probe_records), type="percentage") as bar:
        bar.update(0)
        for count, probe_record in enumerate(probe_records):
            blast_result = run_blast(probe_record, args.database_path.name)
            filtered_results = filter_results(blast_result)
            print_results(filtered_results, probe_record["title"], run_time, args.output_dir)
            bar.update(count)
            # break


if __name__ == "__main__":
    run()
