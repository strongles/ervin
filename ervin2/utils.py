from datatypes import Fasta
from datetime import datetime


def read_file_data(filename):
    with open(filename) as file_in:
        for line in file_in.readlines():
            yield line.strip()


def filename_timestamp():
    current_timestamp = datetime.now()
    return current_timestamp.strftime("%Y-%m-%d_%H-%M-%S")


def generate_output_filepath(base, run_timestamp, extension="fasta"):
    return f"{base.replace('<', '')}_{run_timestamp}.{extension}"


def fasta_file_records(filename):
    title = ""
    seq = []
    for line in read_file_data(filename):
        if line.startswith(">"):
            if title:
                yield Fasta(title, seq)
            title = line
            seq = []
        else:
            seq.append(line)
    if title and seq:
        yield Fasta(title, seq)


def write_tsv_file(records, source_accession_id, filename):
    with open(filename, 'w') as file_out:
        for record in records:
            file_out.write(record.to_tsv_record(source_accession_id))


def write_fasta_file(fasta_records: list[Fasta], filename):
    with open(filename, 'w') as file_out:
        for fasta in fasta_records:
            file_out.write(str(fasta))
