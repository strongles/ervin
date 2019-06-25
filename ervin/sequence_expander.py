from .defaults import CONFIG_FILEPATH
from .defaults import NEWLINE
from .defaults import TEMP_BLASTN_OUTPUT
from .defaults import TEMP_FASTA_FILE
from .defaults import LTR_LOWER
from .defaults import LTR_UPPER
from .defaults import LTR_OUTFILE
from .defaults import ENV_UPPER

from .exceptions import BadConfigFormatException, IncompleteArgsException

from .scaf_file import ScafRecord

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

import argparse
import json
import os


REQUIRED_ARGS = ["file"]
#     [
#     "start",
#     "end",
#     "accession_id",
#     "range_expansion_size",
#     "range_size"
# ]


def parse_input_file(filepath):
    if os.path.isfile(filepath):
        file_entries = []
        with open(filepath) as input_file:
            for line in input_file:
                line_tokens = line.strip().split(" ")
                file_entries.append({
                    "accession_id": line_tokens[0],
                    "start": int(line_tokens[1]),
                    "end": int(line_tokens[2])
                })
        return file_entries
    else:
        raise FileExistsError(f"File {filepath} does not exist.")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        help="File containing genome segment definitions",
                        type=str,
                        required=True),
    parser.add_argument("-res", "--range_expansion_size",
                        help="the number of positions by which to expand "
                             "the selection window from the genome",
                        type=int)
    parser.add_argument("-ff", "--fasta_filepath",
                        help="filepath to the fasta db to query for the full genome sequence",
                        type=str)
    return parser.parse_args()
    # parser = argparse.ArgumentParser()
    # parser.add_argument("-s", "--start",
    #                     help="start of the range to examine in the genome",
    #                     type=int,
    #                     required=True)
    # parser.add_argument("-e", "--end",
    #                     help="end of the range to examine in the genome",
    #                     type=int,
    #                     required=True)
    # parser.add_argument("-aid", "--accession_id",
    #                     help="accession id",
    #                     type=str,
    #                     required=True)
    # parser.add_argument("-res", "--range_expansion_size",
    #                     help="the number of positions by which to expand the
    #                     selection window from the genome",
    #                     type=int)
    # parser.add_argument("-rs", "--range_size",
    #                     help="TODO: Figure out what this is for",
    #                     type=int)
    # parser.add_argument("-f", "--fasta_filepath",
    #                     help="filepath to the fasta db to query for the full genome sequence",
    #                     type=str)
    # return parser.parse_args()


def get_conf_settings():
    if os.path.isfile(CONFIG_FILEPATH):
        with open(CONFIG_FILEPATH) as config_file:
            try:
                return json.load(config_file)
            except Exception:
                raise BadConfigFormatException(Exception)
    else:
        print("Config file not present. All required arguments must be passed manually")
        return {}


def validate_settings(conf_settings, parsed_args):
    validation = {arg: (arg not in conf_settings and arg not in parsed_args)
                  for arg in REQUIRED_ARGS}
    print(validation)
    if any([(arg not in conf_settings and arg not in parsed_args) for arg in REQUIRED_ARGS]):
        raise IncompleteArgsException(f"Args missing. Ensure all are provided: "
                                      f"\n{NEWLINE.join(REQUIRED_ARGS)}")
    return {**conf_settings, **vars(parsed_args)}


def compute_ranges(range_start, range_end):
    return range_start - 5000, range_end + 5000


def get_from_db(range_start, range_end, accession_id, filepath):
    # TODO: Find out DB interface
    for seq_record in SeqIO.parse(filepath, "fasta"):
        if seq_record.id == accession_id:
            return seq_record.seq[range_start:range_end]
    raise Exception("Accession ID not found.")


def run_blast_against_tempfile():
    # outfmt=15 -> JSON
    # outfmt=5 -> XML
    blast_command = NcbiblastnCommandline(out=TEMP_BLASTN_OUTPUT,
                                          query=TEMP_FASTA_FILE,
                                          subject=TEMP_FASTA_FILE,
                                          outfmt=15,
                                          strand="both")
    blast_command()


def write_result_to_tempfile(record):
    record.print_to_file(TEMP_FASTA_FILE)


def read_blast_result_from_file():
    with open(TEMP_BLASTN_OUTPUT) as blast_file:
        return json.load(blast_file)


def get_blast_results():
    blast_data = read_blast_result_from_file()["BlastOutput2"][0]["report"]
    return blast_data["results"]["bl2seq"][0]


def parse_blast_fasta_title(title):
    return_dict = {}
    tokens = title.split(" ")
    return_dict["acc_id"] = tokens[0]
    return_dict["start"] = int(tokens[1])
    return_dict["end"] = int(tokens[2])
    return return_dict


def determine_ltr_hits(blast_results, db_path):
    # print("In determine LTRs")
    ltr_return = []
    blast_score_groups = {}
    for hit in blast_results["hits"][0]["hsps"]:
        if hit["bit_score"] not in blast_score_groups:
            blast_score_groups[hit["bit_score"]] = [hit]
        else:
            blast_score_groups[hit["bit_score"]].append(hit)
        # print(hit)
    # print(f"Blast score groups: {blast_score_groups}")
    potential_ltrs = []
    for key, value in blast_score_groups.items():
        print(f"score: {key} hits: {len(value)}")
        if len(value) > 1:
            if LTR_LOWER < value[0]["align_len"] < LTR_UPPER:
                potential_ltrs.append(value)
    # print(f"Potential LTRS: {potential_ltrs}")
    for potential_ltr in potential_ltrs:
        start = 0
        end = 0
        for record in potential_ltr:
            start = min(record["query_from"], record["hit_from"])
            end = max(record["query_to"], record["hit_to"])
            print(f"Sequence from {start} to {end}; length: {end-start}")
        parsed_title = parse_blast_fasta_title(blast_results["query_title"])
        full_segment_start = start + parsed_title["start"]
        full_segment_end = end + parsed_title["start"]
        # print(f"Original start: {parsed_title[1]} Original end: {parsed_title[2]}")
        # print(f"Internal start: {start_point} Internal end: {end_point}")
        print(f"Extracted sequence length: {start}-{end} = {end- start}")
        print(f"Full sequence start and end: {full_segment_start}-{full_segment_end} "
              f"Sequence length should be {full_segment_end-full_segment_start}")
        if full_segment_end - full_segment_start <= ENV_UPPER:
            scaf_record = ScafRecord(accession_id=parsed_title["acc_id"],
                                     first_position=full_segment_start,
                                     second_position=full_segment_end,
                                     segment=get_from_db(full_segment_start,
                                                         full_segment_end,
                                                         parsed_title["acc_id"],
                                                         db_path))
            ltr_return.append(str(scaf_record))
    return ltr_return


def print_ltr_candidates_to_file(ltr_list):
    with open(LTR_OUTFILE, "w") as outfile:
        for ltr in ltr_list:
            outfile.write(ltr)


def run():
    args = validate_settings(get_conf_settings(), parse_args())
    # acc_id = args["accession_id"]
    # start = args["start"]
    # end = args["end"]
    records = parse_input_file(args["file"])
    db_filepath = args["fasta_filepath"]
    ltr_candidates = []
    for record in records:
        new_range_start, new_range_end = compute_ranges(record["start"], record["end"])
        scaf_record = ScafRecord(accession_id=record["accession_id"],
                                 first_position=new_range_start,
                                 second_position=new_range_end,
                                 segment=get_from_db(new_range_start, new_range_end,
                                                     record["accession_id"], db_filepath))
        write_result_to_tempfile(scaf_record)
        run_blast_against_tempfile()
        blast_results = get_blast_results()
        ltr_candidates.extend(determine_ltr_hits(blast_results, db_filepath))
    print_ltr_candidates_to_file(ltr_candidates)
    # TODO: Plug the LTR segments into the VIRUS DB to confirm them further


if __name__ == "__main__":
    run()
