import json
from Bio.Blast.Applications import NcbiblastnCommandline
from config import get_config


def _ncbi_blast_commandline_wrapper(**kwargs):
    command = NcbiblastnCommandline(**kwargs)
    # Command execution returns a tuple, the actual result is the first element
    command_result = json.loads(command()[0])
    # We don't care about the first couple of layers of the JSON,
    # so we can pare it down to the interesting stuff here
    return command_result["BlastOutput2"][0]["report"]["results"]["search"]


def _tblastn_e_value_limited(query_file, db_path, e_value):
    return _ncbi_blast_commandline_wrapper(cmd="tblastn",
                                           outfmt=15,
                                           query=query_file,
                                           db=db_path,
                                           evalue=e_value)


def _tblastn_simple(query_file, db_path):
    return _ncbi_blast_commandline_wrapper(cmd="tblastn",
                                           outfmt=15,
                                           query=query_file,
                                           db=db_path)


def tblastn(query_file, db_name, e_value=None):
    config = get_config()
    db_path = f"{config.genome_db_storage}{db_name}"
    if e_value:
        return _tblastn_e_value_limited(query_file, db_path, e_value)
    else:
        return _tblastn_simple(query_file, db_path)
