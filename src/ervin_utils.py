from pathlib import Path
import datetime

DEFAULT_OUTPUT_DIR = Path.cwd() / "OUTPUT"
TEMP_FASTA_FILE = Path("/tmp/temp.fasta")
TEMP_TBLASTN_OUTPUT = Path("/tmp/temp_tblastn.tsv")
NEWLINE = "\n"


def format_timestamp_for_filename():
    current_timestamp = datetime.datetime.now()
    return current_timestamp.strftime("%Y-%m-%d_%H-%M-%S")


def read_and_sanitise_raw_data(filename):
    with open(filename) as fasta_input:
        raw_data = fasta_input.readlines()
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


def print_to_fasta_file(filename, fasta_list):
    with open(filename, 'w') as fasta_out:
        for fasta_record in fasta_list:
            fasta_out.write(f"{fasta_record['title']}{NEWLINE}{fasta_record['seq']}{NEWLINE}")
