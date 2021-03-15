import logging
import json
import os
from tqdm import tqdm

from blast_wrapper import tblastn
from utils import filename_timestamp, write_fasta_file, write_tsv_file
from utils import fasta_file_records
from utils import generate_output_filepath
from datatypes import BlastResult
from datatypes import BlastHit
from datatypes import Fasta
from config import get_config


ALIGNMENT_THRESHOLD_LEN = None
LOGGER = logging.getLogger("ervin/blast_probes")
LOGGER.setLevel(logging.INFO)


def default_output_dir():
    return os.getcwd()


def temp_fasta_file(fasta_record: Fasta):
    temp_filepath = "/tmp/temp.fasta"
    with open(temp_filepath, "w") as fasta_out:
        LOGGER.info(f"Writing temporary FASTA file to {temp_filepath}")
        fasta_out.write(str(fasta_record))
    return temp_filepath


def run_blast(probe_fasta: Fasta, genome_db, e_val_limit):
    LOGGER.info(f"Running BLAST with args: {genome_db=}, {e_val_limit=}")
    blast_result = tblastn(temp_fasta_file(probe_fasta), genome_db, e_val_limit)
    return BlastResult(**blast_result)


def _exceeds_length_threshold(blast_hit: BlastHit):
    return blast_hit.alignment_length > ALIGNMENT_THRESHOLD_LEN


def filtered_blast_hits(blast_hits: list[BlastHit]):
    filters = [
        _exceeds_length_threshold
    ]
    for blast_hit in blast_hits:
        if all(filter_case(blast_hit) for filter_case in filters):
            yield blast_hit


def run_blast_probes(filename, genome_db, output_dir, align_thresh, e_val_limit, run_ts=None):
    config = get_config()
    global ALIGNMENT_THRESHOLD_LEN
    ALIGNMENT_THRESHOLD_LEN = align_thresh if align_thresh else config.probe_blaster_align_len_thresh

    file_ts = run_ts if run_ts else filename_timestamp()
    LOGGER.info(f"Run time: {run_ts}")
    fasta_records = fasta_file_records(filename)
    output_filepaths = []

    for fasta in tqdm(fasta_records):
        blast_result = run_blast(fasta, genome_db, e_val_limit)
        filtered_hits = filtered_blast_hits(blast_result.hit_results)
        output_filepath = generate_output_filepath(
            os.path.join(output_dir, fasta.title.replace(">", "")),
            file_ts,
            extension="tsv"
        )
        write_tsv_file(filtered_hits, blast_result.query_accession_id, output_filepath)
        output_filepaths.append(output_filepath)
