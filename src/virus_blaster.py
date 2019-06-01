from pathlib import Path

from ervin_utils import DEFAULT_OUTPUT_DIR, read_from_fasta_file, \
    format_timestamp_for_filename, print_to_fasta_file
from Bio.Blast import NCBIWWW, NCBIXML
import progressbar
import argparse


VIRUS_DATABASE_ID = "10239"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir",
                        help="Directory in which to create the output files",
                        type=str,
                        required=False)
    parser.add_argument("-f", "--filelist",
                        help="List of files to run against the virus database",
                        type=argparse.FileType('r'),
                        nargs='+',
                        required=True)
    return parser.parse_args()


def get_data_from_all_files(args):
    file_data = []
    for filepath in args.filelist:
        file_data.append(read_from_fasta_file(filepath.name))
    return file_data


def get_blast_results(data):
    result_set = {}
    with progressbar.ProgressBar(max_value=len(data), type="percentage") as bar:
        for enum, datum in enumerate(data):
            bar.update(enum)
            result_handle = NCBIWWW.qblast("tblastn", "refseq_genomic",
                                           datum["seq"], hitlist_size=1,
                                           entrez_query=f"txid{VIRUS_DATABASE_ID}[ORGN]")
            # Due to the hitlist size sent to the API,
            # we only care about the very first result in this
            blast_result = next(NCBIXML.parse(result_handle))
            hit_descriptions = blast_result.descriptions
            if len(hit_descriptions) > 0:
                result_title = hit_descriptions[0].title
            else:
                result_title = "not_found"
            if result_title not in result_set:
                result_set[result_title] = [datum]
            else:
                result_set[result_title].append(datum)

        return result_set


def print_virus_results_to_file(virus_name, result_list, run_ts):
    filename = Path(DEFAULT_OUTPUT_DIR) / f"{virus_name}_{run_ts}.fasta"
    print_to_fasta_file(filename, result_list)


def run():
    args = parse_args()
    run_stamp = format_timestamp_for_filename()
    input_data = get_data_from_all_files(args)
    with progressbar.ProgressBar(max_value=len(input_data), type="percentage") as outer_bar:
        for enum, file_data in enumerate(input_data):
            outer_bar.update(enum)
            blast_results = get_blast_results(file_data)
            for result_set in blast_results:
                print_virus_results_to_file(result_set, blast_results[result_set], run_stamp)


if __name__ == "__main__":
    run()
