class ScafRecord:
    """Class for holding a single line-record as represented in a scaf.txt file"""
    accession_id: str
    start: int
    end: int
    direction: str
    segment: str

    def __init__(self,
                 record_line=None,
                 accession_id=None,
                 first_position=None,
                 second_position=None,
                 segment=None):
        if record_line is not None:
            line_tokens = record_line.strip().split("\t")
            accession_id = line_tokens[0]
            first_position = int(line_tokens[1])
            second_position = int(line_tokens[2])
            segment = line_tokens[3]
        self.accession_id = accession_id
        if first_position < second_position:
            self.start = first_position
            self.end = second_position
            self.direction = "P"
            self.segment = segment
        else:
            self.start = second_position
            self.end = first_position
            self.direction = "N"
            self.segment = segment[::-1]

    def __str__(self):
        return f">{self.accession_id} {self.start} {self.end} {self.direction}\n{self.segment}\n"

    def print_to_file(self, filepath):
        with open(filepath, "w") as outfile:
            outfile.write(str(self))
