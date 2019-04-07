from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from src.exceptions import *
from src.defaults import *
from src.scaf_file import ScafRecord
import argparse
import json
import os


REQUIRED_ARGS = [
    "start",
    "end",
    "accession_id",
    "range_expansion_size",
    "range_size"
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--start",
                        help="start of the range to examine in the genome",
                        type=int,
                        required=True)
    parser.add_argument("-e", "--end",
                        help="end of the range to examine in the genome",
                        type=int,
                        required=True)
    parser.add_argument("-aid", "--accession_id",
                        help="accession id",
                        type=str,
                        required=True)
    parser.add_argument("-res", "--range_expansion_size",
                        help="the number of positions by which to expand the selection window from the genome",
                        type=int)
    parser.add_argument("-rs", "--range_size",
                        help="TODO: Figure out what this is for",
                        type=int)
    return parser.parse_args()


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
    validation = {arg: (arg not in conf_settings and arg not in parsed_args) for arg in REQUIRED_ARGS}
    print(validation)
    if any([(arg not in conf_settings and arg not in parsed_args) for arg in REQUIRED_ARGS]):
        raise IncompleteArgsException(f"Args missing. Ensure all are provided: \n{NEWLINE.join(REQUIRED_ARGS)}")
    # parsed_args.__dict__()
    return {**conf_settings, **vars(parsed_args)}
    # return dict(conf_settings, **parsed_args)


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
        return json.load(blast_file.read())


def get_blast_results():
    # Filering out some of the bumf of unecessary nesting here, as nothing before the report data is really meaningful
    # I'm looking at you XML
    blast_data = read_blast_result_from_file()["BlastOutput2"]["report"]
    return blast_data["results"]["bl2seq"]


def parse_blast_fasta_title(title):
    return_dict = {}
    tokens = title.split(" ")
    return_dict["acc_id"] = tokens[0]
    return_dict["start"] = tokens[1]
    return_dict["end"] = tokens[2]
    return return_dict


def determine_ltr_hits(blast_results, db_path):
    ltr_return = []
    blast_score_groups = {}
    for hit in blast_results:
        if hit["score"] not in blast_score_groups:
            blast_score_groups[hit["bit_score"]] = [hit]
        else:
            blast_score_groups[hit["bit_score"]].append(hit)
    potential_ltrs = []
    for key, value in blast_score_groups.items():
        if len(value) == 2:
            potential_ltrs.append(value)
    for potential_ltr in potential_ltrs:
        if potential_ltr[0]["query_from"] == potential_ltr[1]["hit_from"] \
                or potential_ltr[0]["hit_from"] == potential_ltr[1]["query_from"]:
            start_pos = min(potential_ltr[0]["query_from"],
                            potential_ltr[0]["hit_from"],
                            potential_ltr[1]["query_from"],
                            potential_ltr[1]["hit_from"])
            end_pos = max(potential_ltr[0]["query_to"],
                          potential_ltr[0]["hit_to"],
                          potential_ltr[1]["query_to"],
                          potential_ltr[1]["hit_to"])
            parsed_title = parse_blast_fasta_title(blast_results["query_title"])
            full_segment_start = start_pos + parsed_title["start"]
            full_segment_end = end_pos + parsed_title["end"]
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
    acc_id = args["accession_id"]
    start = args["start"]
    end = args["end"]
    db_filepath = args["fasta_filepath"]

    new_range_start, new_range_end = compute_ranges(start, end)
    scaf_record = ScafRecord(accession_id=acc_id,
                             first_position=new_range_start,
                             second_position=new_range_end,
                             segment=get_from_db(new_range_start, new_range_end, acc_id, db_filepath))
    write_result_to_tempfile(scaf_record)
    run_blast_against_tempfile()
    blast_results = get_blast_results()
    ltr_candidates = determine_ltr_hits(blast_results, db_filepath)
    print_ltr_candidates_to_file(ltr_candidates)
    # TODO: Plug the LTR segments into the VIRUS DB to confirm them further


if __name__ == "__main__":
    run()
