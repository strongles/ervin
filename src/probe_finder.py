from math import ceil
import os

FASTA_OUTPUT = "output/result.fasta"
TSV_OUTPUT = "output/result.tsv"

INPUT_FILE_ONE = "/Users/rob/PycharmProjects/ervin/data/first_probe_file.tsv"
INPUT_FILE_TWO = "/Users/rob/PycharmProjects/ervin/data/second_probe_file.tsv"


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
            self.alignment_length = line_tokens[6]
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
        return self.scaffold < other.scaffold

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
                               self.e_value, self.alignment_length, self.acc_sequence, self.scaffold_alignment,
                               self.frame]]
        tab = "\t"
        return f"{tab.join(stringified_output)}\n"

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
                                                       "N" * ceil((second.start - first.end) / 3),
                                                       second.scaffold_alignment])
        elif first.end > second.start:
            overlap = ceil((first.end - second.start) / 3)  # Account for nucleotide -> protein conversion
            overrides["scaffold_alignment"] = "".join([first.scaffold_alignment,
                                                       second.scaffold_alignment[:overlap]])
        align_length = len(overrides["scaffold_alignment"]) - str.count("-", overrides["scaffold_alignment"])
        overrides["alignment_length"] = align_length
        overrides["accession_id"] = first.accession_id + "_" + second.accession_id
        return ProbeData(first, overrides)


def read_probe_records_from_file(filename):
    record_list = []
    with open(filename) as probe_in:
        for line in probe_in:
            record_list.append(ProbeData(line))
    record_dict = {}
    for record in record_list:
        if record.scaffold not in record_dict:
            record_dict[record.scaffold] = [record]
        else:
            record_dict[record.scaffold].append(record)
    return record_dict


# Account for incomplete probe, merge the alignment so that the
def is_near_neighbour(master, comparitor):
    first_diff = master.start - comparitor.end
    second_diff = comparitor.start - master.end
    return (0 < first_diff <= 50 or 0 < second_diff <= 50) and master.direction == comparitor.direction and master.frame == comparitor.frame


def is_range_extension(master, comparitor):
    return (((master.start <= comparitor.start) and (comparitor.start <= master.end < comparitor.end)) or
            ((comparitor.start <= master.start) and (master.start <= comparitor.end < master.end)) or
            ((master.end >= comparitor.end) and (comparitor.start <= master.start < comparitor.end)) or
            ((comparitor.end >= master.end) and (master.start <=comparitor.start < master.end))) and \
            master.direction == comparitor.direction and master.frame == comparitor.frame


def is_superset(master, comparitor):
    return master.start >= comparitor.start and master.end <= comparitor.end and master.frame == comparitor.frame


def write_to_files(record):
    with open(FASTA_OUTPUT, "a") as fasta:
        with open(TSV_OUTPUT, "a") as tsv:
            fasta.write(record.to_fasta())
            tsv.write(record.to_tsv())


def set_up_output_files():
    with open(FASTA_OUTPUT, "w"):
        with open(TSV_OUTPUT, "w"):
            pass


def unique_scaffolds(source_one, source_two):
    source_one_uniques = [scaffold for scaffold in source_one.keys() if scaffold not in source_two.keys()]
    source_two_uniques = [scaffold for scaffold in source_two.keys() if scaffold not in source_one.keys()]
    return_set = set()
    for unique in source_one_uniques:
        return_set.update(source_one[unique])
        del source_one[unique]
    for unique in source_two_uniques:
        return_set.update(source_two[unique])
        del source_two[unique]
    return return_set


def run():
    set_up_output_files()
    first_probe_data = read_probe_records_from_file(INPUT_FILE_ONE)
    second_probe_data = read_probe_records_from_file(INPUT_FILE_TWO)
    # Add any scaffolds which don't appear in the other file automatically
    output_set = unique_scaffolds(first_probe_data, second_probe_data)

    for scaffold, records in first_probe_data.items():
        for record in records:
            current_print_candidate = record
            for comparitor in second_probe_data[scaffold]:
                if is_superset(current_print_candidate, comparitor):
                    current_print_candidate = comparitor
                elif is_near_neighbour(current_print_candidate, comparitor) \
                        or is_range_extension(current_print_candidate, comparitor):
                    current_print_candidate = ProbeData.merge_records(current_print_candidate, comparitor)
            output_set.add(current_print_candidate)

    [write_to_files(output) for output in sorted(list(output_set))]


if __name__ == "__main__":
    run()
