"""
for row in file1:
    for line in file2:
        if accession_ids match and scaffolds match:

"""

FASTA_OUTPUT = "output/result.fasta"
TSV_OUTPUT = "output/result.tsv"


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
    unknown_data = str
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
            if overrides is not None:
                for attribute, value in overrides:
                    setattr(self, attribute, value)

        elif source is not None:
            line_tokens = source.strip().split("\t")
            self.accession_id = line_tokens[0]
            self.scaffold = line_tokens[1]
            self.scaffold_length = int(line_tokens[2])
            self.e_value = line_tokens[5]
            self.alignment_length = line_tokens[6]
            self.acc_sequence = line_tokens[7]
            self.scaffold_alignment = line_tokens[8]
            self.frame = line_tokens[9]
            self.matched = False
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

    def to_fasta(self):
        return f">{self.scaffold} {self.start} {self.end} {self.direction}\n{self.scaffold_alignment}\n"

    def to_tsv(self):
        if self.direction == "P":
            start = self.start
            end = self.end
        else:
            start = self.end
            end = self.start
        stringified_output = [str(output) for output in
                              [self.accession_id, self.scaffold, self.scaffold_length, start, end,
                               self.e_value, self.alignment_length, self.acc_sequence, self.scaffold_alignment,
                               self.frame]]
        tab = "\t"
        return f"{tab.join(stringified_output)}\n"

    @staticmethod
    def merge_records(a, b):
        if a.end > b.start:
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
                                                       "N" * (second.start - first.end / 3),
                                                       second.scaffold_alignment])
        elif first.end > second.start:
            overlap = first.end - second.start / 3  # Account for nucleotide -> protein conversion
            overrides["scaffold_alignment"] = "".join([first.scaffold_alignment,
                                                       second.scaffold_alignment[:overlap]])
        align_length = len(overrides["scaffold_alignment"]) - str.count("-", overrides["scaffold_alignment"])
        overrides["alignment_length"] = align_length
        return ProbeData(first, overrides)


def read_probe_records_from_file(filename):
    return_list = []
    with open(filename) as probe_in:
        for line in probe_in:
            return_list.append(ProbeData(line))
    return return_list


# Account for incomplete probe, merge the alignment so that the
def is_near_neighbour(master, comparitor):
    return (master.start - comparitor.end <= 50 or comparitor.start - master.end <= 50) \
           and master.direction == comparitor.direction


def is_range_extension(master, comparitor):
    return (comparitor.start < master.start and master.end < comparitor.end) or \
           (master.start > comparitor.start and master.end > comparitor.end)


def is_superset(master, comparitor):
    return master.start > comparitor.start and master.end < comparitor.end


def write_to_files(record):
    with open(FASTA_OUTPUT, "a") as fasta:
        with open(TSV_OUTPUT, "a") as tsv:
            fasta.write(record.to_fasta())
            tsv.write(record.to_tsv())


def clear_output_files():
    with open(FASTA_OUTPUT, "a"):
        with open(TSV_OUTPUT, "a"):
            pass


def run():
    clear_output_files()
    first_probe_data = read_probe_records_from_file("first_probe_file")
    second_probe_data = read_probe_records_from_file("second_probe_file")
    for probe_record in first_probe_data:
        for probe_comparitor in second_probe_data:
            if probe_record.scaffold == probe_comparitor.scaffold:
                probe_record.matched = True
                probe_comparitor.matched = True
                if is_superset(probe_record, probe_comparitor):
                    write_to_files(probe_comparitor)
                elif is_near_neighbour(probe_record, probe_comparitor) \
                        or is_range_extension(probe_record, probe_comparitor):
                    write_to_files(ProbeData.merge_records(probe_record, probe_comparitor))
                else:
                    continue
        if not probe_record.matched:
            write_to_files(probe_record)
    unmatched_comparitors = [comparitor for comparitor in second_probe_data if comparitor.matched is False]
    map(write_to_files, unmatched_comparitors)


if __name__ == "__main__":
    run()
