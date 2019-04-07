from Bio import SeqIO
import argparse
import json
import os
# from Bio.Blast.NCBIWWW import qblast


class BadConfigFormatException(Exception):
    pass


class IncompleteArgsException(Exception):
    pass


FILEPATH = "../data/GCA_004026985.1_MyoMyo_v1_BIUU_genomic.fna"
CONFIG_FILEPATH = "config/config.json"
NEWLINE = "\n"
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


if __name__ == "__main__":
    args = validate_settings(get_conf_settings(), parse_args())
    acc_id = args["accession_id"]
    start = args["start"]
    end = args["end"]
    db_filepath = args["fasta_filepath"]

    new_range_start, new_range_end = compute_ranges(start, end)
    result = get_from_db(new_range_start, new_range_end, acc_id, db_filepath)
