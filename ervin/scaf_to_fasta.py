from .scaf_file import ScafRecord

input_filepath = "scaf.txt"
output_filepath = "fasta_output.txt"
fasta_filepath = "blah.fasta"


def construct_file_line(scaf_record):
    return f">{scaf_record.accession_id}_{scaf_record.start}_{scaf_record.end}_" \
        f"{scaf_record.direction}\n{scaf_record.segment}\n"


if __name__ == "__main__":
    with open(input_filepath) as infile:
        with open(output_filepath, "w") as outfile:
            for line in infile:
                scaf_record = ScafRecord(record_line=line)
                outfile.write(construct_file_line(scaf_record))
