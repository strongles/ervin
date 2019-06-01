from ervin_utils import format_timestamp_for_filename, DEFAULT_OUTPUT_DIR
from exceptions import InvalidPathException
from probe_data import ProbeData
import progressbar
import argparse
import os


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


def write_to_files(output_files, record):
    (fasta_output, tsv_output) = output_files
    with open(fasta_output, "a") as fasta:
        with open(tsv_output, "a") as tsv:
            fasta.write(record.to_fasta())
            tsv.write(record.to_tsv())


def set_up_output_files(output_dir):
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
    source_one_uniques = [
        scaffold for scaffold in source_one.keys() if scaffold not in source_two.keys()
    ]
    source_two_uniques = [
        scaffold for scaffold in source_two.keys() if scaffold not in source_one.keys()
    ]
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


def is_align_length_threshold_satisfied(args, record):
    if args.alignment_len_threshold is None:
        return True
    elif record.alignment_length > args.alignment_len_threshold:
        return True
    else:
        return False


def find_probes(args, first_probe_data, second_probe_data):
    output_data = unique_scaffolds(first_probe_data, second_probe_data)

    for scaffold, records in first_probe_data.items():
        for record in records:
            if is_align_length_threshold_satisfied(args, record):
                current_print_candidate = record
                for comparitor in second_probe_data[scaffold]:
                    if current_print_candidate.is_superset(comparitor):
                        current_print_candidate = comparitor
                    elif current_print_candidate.is_near_neighbour(comparitor) \
                            or current_print_candidate.is_range_extension(comparitor):
                        current_print_candidate = ProbeData.merge_records(current_print_candidate,
                                                                          comparitor)
                if scaffold not in output_data:
                    output_data[scaffold] = {current_print_candidate}
                else:
                    output_data[scaffold].add(current_print_candidate)
    return output_data


def track_recursive_probe_finder_progress(func):
    try:
        print(len(progressbars))
    except NameError:
        progressbars = {}

    def wrapper(file_list, **kwargs):
        try:
            progressbars[func].update(progressbars[func].max_value - len(file_list))
        except KeyError:
            progressbars[func] = progressbar.ProgressBar(max_value=len(file_list),
                                                         type="percentage",
                                                         prefix="Total Run: ")
            progressbars[func].update(0)
        return func(file_list, **kwargs)
    return wrapper


@track_recursive_probe_finder_progress
def find_probes_recursively(file_list, args, tail=None):
    if tail is None:
        first_probe_data = read_probe_records_from_file(file_list[0])
        second_probe_data = read_probe_records_from_file(file_list[1])

        if len(file_list) == 2:
            return find_probes(args, first_probe_data, second_probe_data)
        elif len(file_list) < 2:
            return find_probes_recursively(file_list[2:],
                                           args,
                                           tail=find_probes(args,
                                                            first_probe_data,
                                                            second_probe_data))
    else:
        probe_data = read_probe_records_from_file(file_list[0])

        if len(file_list) == 1:
            return find_probes(args, tail, probe_data)
        else:
            return find_probes_recursively(file_list[1:], args,
                                           tail=find_probes(args, tail, probe_data))


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
            result = find_probes_recursively(input_files, args)
    else:
        file_list = read_filenames_from_manifest(args.manifest.name)
        result = find_probes_recursively(file_list, args)

    if result is not None:
        result_data = flatten_results(result)
        output_files = set_up_output_files(args.output_dir)
        [write_to_files(output_files, output) for output in result_data]
    else:
        raise Exception("No results after running probe_finder.")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir",
                        help="Directory in which to create the output files",
                        type=str,
                        default=DEFAULT_OUTPUT_DIR,
                        required=False)
    parser.add_argument("-a", "--alignment_len_threshold",
                        help="Minimum length threshold that BLAST result "
                             "alignment sequence lengths should exceed",
                        type=int,
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
