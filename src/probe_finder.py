from exceptions import InvalidPathException, IncompleteArgsException
from math import ceil
import argparse
import os
import datetime


DEFAULT_OUTPUT_DIR = os.path.join(os.getcwd(), "OUTPUT")


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


def write_to_files(output_files, record):
    (fasta_output, tsv_output) = output_files
    with open(fasta_output, "a") as fasta:
        with open(tsv_output, "a") as tsv:
            fasta.write(record.to_fasta())
            tsv.write(record.to_tsv())


def format_timestamp_for_filename():
    current_timestamp = datetime.datetime.now()
    return current_timestamp.strftime("%Y-%m-%d_%H:%M:%S")


def set_up_output_files(output_dir):
    if output_dir is None:
        output_dir = DEFAULT_OUTPUT_DIR
    if not os.path.isdir(output_dir):
        if os.path.exists(output_dir):
            raise InvalidPathException(f"Invalid output path provided: {output_dir}")
        else:
            os.makedirs(output_dir)
    run_time = format_timestamp_for_filename()
    fasta_output_path = os.path.join(output_dir, f"probe_finder-{run_time}.fasta")
    tsv_output_path = os.path.join(output_dir, f"probe_finder-{run_time}.tsv")
    with open(fasta_output_path, "w"):
        with open(tsv_output_path, "w"):
            pass
    return fasta_output_path, tsv_output_path


def unique_scaffolds(source_one, source_two):
    source_one_uniques = [scaffold for scaffold in source_one.keys() if scaffold not in source_two.keys()]
    source_two_uniques = [scaffold for scaffold in source_two.keys() if scaffold not in source_one.keys()]
    return_dict = {}
    for unique in source_one_uniques:
        return_dict[unique] = source_one[unique]
        del source_one[unique]
    for unique in source_two_uniques:
        return_dict[unique] = source_two[unique]
        del source_two[unique]
    return return_dict


def flatten_results(result_data):
    result_list = []
    for _, value in result_data.items():
        result_list.extend(list(value))
    return sorted(result_list)


def find_probes(first_probe_data, second_probe_data):
    # first_probe_data = read_probe_records_from_file(input_file_one)
    # second_probe_data = read_probe_records_from_file(input_file_two)
    # Add any scaffolds which don't appear in the other file automatically
    output_data = unique_scaffolds(first_probe_data, second_probe_data)

    for scaffold, records in first_probe_data.items():
        for record in records:
            current_print_candidate = record
            for comparitor in second_probe_data[scaffold]:
                if is_superset(current_print_candidate, comparitor):
                    current_print_candidate = comparitor
                elif is_near_neighbour(current_print_candidate, comparitor) \
                        or is_range_extension(current_print_candidate, comparitor):
                    current_print_candidate = ProbeData.merge_records(current_print_candidate, comparitor)
            if scaffold not in output_data:
                output_data[scaffold] = {current_print_candidate}
            else:
                output_data[scaffold].add(current_print_candidate)
    return output_data


def find_probes_recursively(file_list, tail=None):
    if tail is None:
        first_probe_data = read_probe_records_from_file(file_list[0])
        second_probe_data = read_probe_records_from_file(file_list[1])

        if len(file_list) == 2:
            return find_probes(first_probe_data, second_probe_data)
        elif len(file_list) < 2:
            return find_probes_recursively(file_list[2:], tail=find_probes(first_probe_data, second_probe_data))
    else:
        probe_data = read_probe_records_from_file(file_list[0])

        if len(file_list) == 1:
            return find_probes(tail, probe_data)
        else:
            return find_probes_recursively(file_list[1:], tail=find_probes(tail, probe_data))


def read_filenames_from_manifest(manifest):
    with open(manifest, 'r') as manifest_file:
        return [line.strip() for line in manifest_file.readlines()]


def run():
    args = parse_args()
    result = None
    if args.file_list is not None:
        input_files = [input_file.name for input_file in args.file_list]
        if len(input_files) == 1:
            raise Exception("Uneccessary run with only one file provided.")
        elif len(input_files) > 1:
            result = find_probes_recursively(input_files)
    else:
        file_list = read_filenames_from_manifest(args.manifest.name)
        result = find_probes_recursively(file_list)

    if result is not None:
        result_data = flatten_results(result)
        output_files = set_up_output_files(args.output_dir)
        [write_to_files(output_files, output) for output in result_data]
    else:
        raise Exception("No results after running probe_finder.")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir",
                        help=str,
                        required=False)
    file_sourcing = parser.add_mutually_exclusive_group(required=True)
    file_sourcing.add_argument("-f", "--file_list",
                               help="Input file list",
                               type=argparse.FileType('r'),
                               nargs='*')
    file_sourcing.add_argument("-m", "--manifest",
                               help="Filepath for a manifest file containing a list of input files",
                               type=argparse.FileType('r')),
    return parser.parse_args()


if __name__ == "__main__":
    run()
