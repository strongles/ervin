from .exceptions import ArgNotSupportedError
from .templates.config_template import CONFIG_TEMPLATE

from collections import namedtuple
from functools import lru_cache
from pathlib import Path

import gzip
import json
import logging
import shutil
import datetime

TEMP_PROBE_BLASTER = "/tmp/probe_blaster"
TEMP_PROBE_FINDER = "/tmp/probe_finder"
TEMP_VIRUS_BLASTER = "/tmp/virus_blaster"
DEFAULT_OUTPUT_DIR = Path.cwd() / "OUTPUT"
TEMP_FASTA_FILE = "/tmp/temp.fasta"
TEMP_TBLASTN_OUTPUT = "/tmp/temp_tblastn.tsv"
NEWLINE = "\n"
VIRUS_DB_SERVER = "ftp.ncbi.nlm.nih.gov"
VIRUS_DB_SERVER_DIR = "refseq/release/viral/"
MAKE_BLASTDB_CMD = "makeblastdb -in '{db_files}' -title {db_name} -out {out_path} -dbtype nucl"
# Gives us a handle to the ERViN home directory to access things like config files
# ERVIN_DIR = Path(__file__).parent.parent
CONFIG_PATH = Path.home() / ".ervin/config.json"
REQUIRED_DIRS = [
    "operational_data_storage",
    "virus_db_storage",
    "genome_db_storage"
]
SUPPORTED_CONFIG_ARGS = [
    *REQUIRED_DIRS,
    "range_expansion_size",
    "range_size",
    "probe_blaster_align_len_thresh",
    "probe_blaster_e_value_thresh",
    "log_location"
]

LOGGER = logging.getLogger(Path(__file__).stem)


def decompress_gz_file(filepath):
    dest_file = Path(filepath).with_suffix("")
    LOGGER.info(f"Extracting {filepath} to {dest_file}")
    with gzip.GzipFile(filepath, "rb") as archive:
        with open(dest_file, "wb") as extracted:
            extracted.write(archive.read())
    return str(dest_file)


def decompress_gz_files(file_list):
    return [decompress_gz_file(filepath) for filepath in file_list]


def homify_path(path_string):
    if path_string.startswith("~"):
        return path_string.replace("~", str(Path.home()), 1)
    else:
        return path_string


def ensure_output_dir_exists(output_dir_path):
    if type(output_dir_path) is str:
        destination_dir = Path(output_dir_path)
    else:
        destination_dir = output_dir_path
    if not destination_dir.anchor:
        destination_dir = Path.cwd() / destination_dir
    elif str(destination_dir).startswith("~"):
        destination_dir = Path(homify_path(str(destination_dir)))
    if not destination_dir.exists():
        Path.mkdir(destination_dir, parents=True)
    return destination_dir


def make_config(elem_dict):
    Config = namedtuple("Config", SUPPORTED_CONFIG_ARGS)
    if any([arg not in SUPPORTED_CONFIG_ARGS for arg in elem_dict.keys()]):
        raise ArgNotSupportedError("Config file contains unsupported argument")
    for arg in SUPPORTED_CONFIG_ARGS:
        if arg not in elem_dict:
            elem_dict[arg] = None
    return Config(**elem_dict)


@lru_cache(maxsize=None)
def get_config():
    if not CONFIG_PATH.exists():
        if not CONFIG_PATH.parent.exists():
            CONFIG_PATH.parent.mkdir(parents=True)
        with open(CONFIG_PATH, 'w') as config_writer:
            json.dump(CONFIG_TEMPLATE, config_writer)
    with open(CONFIG_PATH) as config_in:
        config_data = json.load(config_in)
        for directory in REQUIRED_DIRS:
            dir_path = Path(homify_path(config_data[directory]))
            if not dir_path.exists():
                dir_path.mkdir(parents=True)
    return make_config(config_data)


def delete_directory_contents(dir_path):
    for file in dir_path.glob("*"):
        file.unlink()


def sanitise_string(input_string):
    danger_characters = "/?$"
    for character in danger_characters:
        input_string = input_string.replace(character, "_")
    return input_string


def format_timestamp_for_filename():
    current_timestamp = datetime.datetime.now()
    return current_timestamp.strftime("%Y-%m-%d_%H-%M-%S")


def read_and_sanitise_raw_data(filename):
    with open(filename) as file_in:
        raw_data = file_in.readlines()
    stripped_newlines_data = [line.strip() for line in raw_data]
    purged_empty_lines_data = [line for line in stripped_newlines_data if line != ""]
    return purged_empty_lines_data


def read_from_fasta_file(filename):
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


def print_to_fasta_file(filename, fasta_list, mode='w'):
    if mode == 'a':
        if not Path.exists(filename):
            with open(filename, 'w'):
                LOGGER.debug(f"Created output .fasta file: {filename}")
                pass
    with open(filename, mode) as fasta_out:
        for fasta_record in fasta_list:
            fasta_out.write(f"{fasta_record['title']}{NEWLINE}{fasta_record['seq']}{NEWLINE}")


def total_result_records(output_filepaths):
    total_records = 0
    for filepath in output_filepaths:
        with open(filepath) as output_file:
            total_records += len(output_file.readlines())
    return total_records
