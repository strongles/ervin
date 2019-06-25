from Bio import SeqIO

input_filepath = "scaf.txt"
output_filepath = "output.txt"
fasta_filepath = "blah.fasta"


class ScafRecord:
    """Class for holding a single line-record as represented in a scaf.txt file"""
    accession_id: str
    start: int
    end: int
    direction: str

    def __init__(self, record_line=None, accession_id=None,
                 first_position=None, second_position=None):
        if record_line is not None:
            line_tokens = record_line.strip().split("\t")
            accession_id = line_tokens[0]
            first_position = int(line_tokens[1])
            second_position = int(line_tokens[2])
        self.accession_id = accession_id
        if first_position < second_position:
            self.start = first_position
            self.end = second_position
            self.direction = "P"
        else:
            self.start = second_position
            self.end = first_position
            self.direction = "N"


def construct_file_line(scaf_record, segment):
    if scaf_record.direction == "N":
        segment = segment[::-1]
    return f">{scaf_record.accession_id}_{scaf_record.start}_{scaf_record.end}_" \
        f"{scaf_record.direction}\n{segment}\n"


def find_sequence_segments():
    with open(input_filepath) as input_file_reader:
        with open(output_filepath, 'w') as output_file_writer:
            scaf_records = {}
            for line in input_file_reader:
                scaf_record = ScafRecord(record_line=line)
                scaf_records[scaf_record.accession_id] = scaf_record
            for seq_record in SeqIO.parse(fasta_filepath, "fasta"):
                if seq_record.id in scaf_records.keys():
                    scaf_record = scaf_records[seq_record.id]
                    output_file_writer.write(
                        construct_file_line(scaf_record,
                                            str(seq_record.seq[scaf_record.start:scaf_record.end])))


if __name__ == "__main__":
    find_sequence_segments()
