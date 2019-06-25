from math import ceil


class ProbeData:
    accession_id = str
    scaffold = str
    scaffold_length = int
    start = int
    end = int
    e_value = str
    alignment_length = int
    acc_sequence = str
    scaffold_alignment = str
    frame = int

    def __init__(self, source=None, overrides=None):
        if isinstance(source, ProbeData):
            self.accession_id = source.accession_id
            self.scaffold = source.scaffold
            self.scaffold_length = source.scaffold_length
            self.e_value = source.e_value
            self.alignment_length = source.alignment_length
            self.acc_sequence = source.acc_sequence
            self.scaffold_alignment = source.scaffold_alignment
            self.frame = source.frame
            self.start = source.start
            self.end = source.end
            self.direction = source.direction
            self.matched = False
            self.printed = False
            if overrides is not None:
                for attribute, value in overrides.items():
                    setattr(self, attribute, value)

        elif source is not None:
            line_tokens = source.strip().split("\t")
            self.accession_id = line_tokens[0]
            self.scaffold = line_tokens[1]
            self.scaffold_length = int(line_tokens[2])
            self.e_value = line_tokens[5]
            self.alignment_length = int(line_tokens[6])
            self.acc_sequence = line_tokens[7]
            self.scaffold_alignment = line_tokens[8]
            self.frame = int(line_tokens[9])
            self.matched = False
            self.printed = False
            first_position = int(line_tokens[3])
            second_position = int(line_tokens[4])
            if first_position < second_position:
                self.start = first_position
                self.end = second_position
                self.direction = "P"
            else:
                self.start = second_position
                self.end = first_position
                self.direction = "N"

    def __hash__(self):
        return hash(repr(self))

    def __lt__(self, other):
        if self.accession_id < other.accession_id:
            return True
        elif self.scaffold < other.scaffold:
            return True
        elif self.start < other.start:
            return True
        elif self.end < other.end:
            return True
        else:
            return False

    def to_fasta(self):
        if self.direction == "P":
            start = self.start
            end = self.end
        else:
            start = self.end
            end = self.start
        return f">{self.scaffold} {start} {end} {self.direction}\n{self.scaffold_alignment}\n"

    def to_tsv(self):
        if self.direction == "P":
            start = self.start
            end = self.end
        else:
            start = self.end
            end = self.start
        stringified_output = [str(output) for output in
                              [self.accession_id, self.scaffold, self.scaffold_length, start, end,
                               self.e_value, self.alignment_length, self.acc_sequence,
                               self.scaffold_alignment, self.frame]]
        tab = "\t"
        return f"{tab.join(stringified_output)}\n"

    def is_superset(self, comparitor):
        return self.start >= comparitor.start \
               and self.end <= comparitor.end \
               and self.frame == comparitor.frame

    def is_near_neighbour(self, comparitor):
        first_diff = self.start - comparitor.end
        second_diff = comparitor.start - self.end
        return (0 < first_diff <= 50
                or 0 < second_diff <= 50) and \
            self.direction == comparitor.direction and \
            self.frame == comparitor.frame

    def is_range_extension(self, comparitor):
        return (((self.start <= comparitor.start)
                 and (comparitor.start <= self.end < comparitor.end))
                or ((comparitor.start <= self.start)
                    and (self.start <= comparitor.end < self.end))
                or ((self.end >= comparitor.end)
                    and (comparitor.start <= self.start < comparitor.end))
                or ((comparitor.end >= self.end)
                    and (self.start <= comparitor.start < self.end))) \
               and self.direction == comparitor.direction \
               and self.frame == comparitor.frame

    def merge_with(self, source):
        return self.merge_records(self, source)

    @staticmethod
    def merge_records(a, b):
        if a.start < b.start:
            first = a
            second = b
        else:
            first = b
            second = a
        overrides = {
            "start": first.start,
            "end": second.end
        }
        if first.end < second.start:
            overrides["scaffold_alignment"] = "".join([first.scaffold_alignment,
                                                       "-" * ceil((second.start - first.end) / 3),
                                                       second.scaffold_alignment])
        elif first.end > second.start:
            overlap = ceil((first.end - second.start) / 3)
            overrides["scaffold_alignment"] = "".join([first.scaffold_alignment,
                                                       second.scaffold_alignment[overlap:]])
        align_length = (overrides["end"] - overrides["start"]) - \
                       overrides["scaffold_alignment"].count("-")
        overrides["alignment_length"] = align_length
        overrides["accession_id"] = first.accession_id + "_" + second.accession_id
        return ProbeData(first, overrides)
