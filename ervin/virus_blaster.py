from .ervin_utils import DEFAULT_OUTPUT_DIR
from .ervin_utils import MAKE_BLASTDB_CMD
from .ervin_utils import TEMP_FASTA_FILE
from .ervin_utils import VIRUS_DB_SERVER
from .ervin_utils import VIRUS_DB_SERVER_DIR
from .ervin_utils import decompress_gz_files
from .ervin_utils import delete_directory_contents
from .ervin_utils import ensure_output_dir_exists
from .ervin_utils import format_timestamp_for_filename
from .ervin_utils import get_config
from .ervin_utils import homify_path
from .ervin_utils import print_to_fasta_file
from .ervin_utils import read_from_fasta_file
from .ervin_utils import sanitise_string
from .ervin_utils import total_result_records

from .probe_blaster import print_probe_to_temp_fasta_file

from Bio.Blast.Applications import NcbiblastnCommandline
from pathlib import Path

import argparse
import ftputil
import json
import logging
import progressbar
import subprocess
import time


VIRUS_DB_DEFAULT = "Viruses"
VIRUS_DATABASE_ID = "10239"
LOGGER = logging.getLogger(Path(__file__).stem)


def get_local_virus_db_version(storage_path):
    version_file = storage_path / "last_updated"
    if not version_file.exists():
        return 0
    else:
        with open(version_file) as version:
            return int(version.read().strip())


def get_newest_ftp_file_timestamp(ftp_handle):
    file_stats = [ftp_handle.stat(file) for file in ftp_handle.listdir(ftp_handle.curdir)]
    return int(max([file_stat.st_mtime for file_stat in file_stats]))


def get_web_virus_db_version():
    with ftputil.FTPHost(VIRUS_DB_SERVER, "anonymous", "") as virus_serv:
        virus_serv.chdir(VIRUS_DB_SERVER_DIR)
        return get_newest_ftp_file_timestamp(virus_serv)


def virus_db_update_required():
    config = get_config()
    virus_db_store = Path(homify_path(config.virus_db_storage))
    if not virus_db_store.exists():
        virus_db_store.mkdir(parents=True)
        LOGGER.info(f"Created virus db directory: {virus_db_store}")
        return True
    else:
        current_version = get_local_virus_db_version(virus_db_store)
        web_version = get_web_virus_db_version()
        return current_version < web_version


def convert_files_to_db(db_name, source_files):
    config = get_config()
    result = subprocess.run(
        MAKE_BLASTDB_CMD.format(
            db_files=" ".join([str(filepath) for filepath in source_files]),
            db_name=db_name,
            out_path=f"{config.virus_db_storage}{db_name}"),
        capture_output=True, shell=True, check=True)
    LOGGER.debug(result)


def update_virus_db():
    if virus_db_update_required():
        LOGGER.info("Newer virus database version found. Updating...")
        virus_db_files = get_virus_db_files()
        decompressed_virus_files = decompress_gz_files(virus_db_files)
        time.sleep(1)
        convert_files_to_db(VIRUS_DB_DEFAULT, decompressed_virus_files)
    else:
        LOGGER.info("Update not required")
        return


def get_virus_db_files():
    config = get_config()
    local_virus_db_path = Path(homify_path(config.virus_db_storage))
    delete_directory_contents(local_virus_db_path)
    files_to_download = []
    downloaded_files = []
    with ftputil.FTPHost(VIRUS_DB_SERVER, "anonymous", "") as virus_serv:
        virus_serv.chdir(VIRUS_DB_SERVER_DIR)
        for file in virus_serv.listdir(virus_serv.curdir):
            if ".fna" in file:
                files_to_download.append(file)
        with progressbar.ProgressBar(max_value=len(files_to_download),
                                     type="percentage",
                                     prefix="Downloading virus DB source files: ") as bar:
            for count, file in enumerate(files_to_download):
                bar.update(count)
                dest = f"{str(local_virus_db_path)}/{file}"
                LOGGER.info(f"Downloading {file} to {dest}")
                virus_serv.download(file, dest)
                downloaded_files.append(dest)
        file_version = get_newest_ftp_file_timestamp(virus_serv)
    with open(local_virus_db_path / "last_updated", 'w') as version_file:
        version_file.write(str(file_version))
    return downloaded_files


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir",
                        help="Directory in which to create the output files",
                        type=str,
                        required=False,
                        default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("-f", "--file",
                        help="File to run against the virus database",
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument("-vdb", "--virus_database",
                        help="Path to a locally-held copy of the Viruses refseq_genome database",
                        type=str,
                        required=False,
                        default=VIRUS_DB_DEFAULT)
    return parser.parse_args()


def get_data_from_file(filepath):
    try:
        return read_from_fasta_file(filepath.name)
    except AttributeError:
        return read_from_fasta_file(filepath)


def get_top_virus_hit(record, db="Viruses"):
    config = get_config()
    print_probe_to_temp_fasta_file(record)
    command = NcbiblastnCommandline(cmd="tblastn",
                                    outfmt=15,
                                    query=TEMP_FASTA_FILE,
                                    db=f"{config.virus_db_storage}{db}")
    LOGGER.debug(command)
    query_result = command()
    parsed_result = json.loads(query_result[0])["BlastOutput2"][0]

    try:
        return parsed_result["report"]["results"]["search"]["hits"][0]["description"][0]["title"]
    except IndexError:
        return "not_found"


def print_virus_results_to_file(virus_name, result_list, run_ts, output_dir):
    destination_dir = ensure_output_dir_exists(output_dir)
    sanitised_virus_name = sanitise_string(virus_name)
    filename = destination_dir / f"{sanitised_virus_name}_{run_ts}.fasta"
    print_to_fasta_file(filename, result_list, mode='a')
    return filename


def viruses_to_total_found(file_list):
    virus_to_count = {}
    for filepath in file_list:
        row_count = total_result_records([filepath])
        virus_name = str(filepath).replace(".fasta", "").split("/")[-1]
        virus_to_count[virus_name] = row_count // 2
    return virus_to_count


def run_virus_blaster(filename=None, db=None, output_dir=None, run_ts=None):
    run_stamp = run_ts if run_ts else format_timestamp_for_filename()
    update_virus_db()
    input_data = get_data_from_file(filename)
    matched_virus_files = set()
    with progressbar.ProgressBar(max_value=len(input_data),
                                 type="percentage",
                                 prefix="Blasting against Viruses: ") as bar:
        for count, file_record in enumerate(input_data):
            bar.update(count)
            top_hit = get_top_virus_hit(file_record, db)
            matched_virus_files.add(print_virus_results_to_file(top_hit, [file_record],
                                                                run_stamp, output_dir))
    return matched_virus_files, viruses_to_total_found(matched_virus_files)


if __name__ == "__main__":
    args = parse_args()
    run_virus_blaster(args.file.name, args.virus_database, args.output_dir)
