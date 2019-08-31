import json
import logging

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

from pathlib import Path

from .ervin_utils import get_config
from .ervin_utils import homify_path

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)

CONFIG = get_config()


class LoggedException(Exception):
    def __init__(self, message):
        LOGGER.error(message)
        super().__init__(message)


class GenomeDbNotFoundException(Exception):
    def __init__(self, message):
        super().__init__(message)


class AccessionIdNotFoundException(Exception):
    def __init__(self, message):
        super().__init__(message)


class Sequence:
    def __init__(self, start, end, sequence, accession_id, genome):
        if start < end:
            self.start = start
            self.end = end
            self.sequence = sequence
            self.direction = "P"
        else:
            self.start = end
            self.end = start
            self.sequence = sequence[::-1]
            self.direction = "N"
        self.accession_id = accession_id
        self.genome = genome


def run(sequence_list, genomes):
    for sequence in sequence_list:
        expanded_sequence = expand_sequence(sequence)
        blast_result = blast_sequence_against_itself(expanded_sequence)
        regions_to_remove = find_regions_to_remove(blast_result)
        removable_region_boundaries = [region for region in map(deletable_region, regions_to_remove)]
        purged_sequence = remove_regions_from_sequence(expanded_sequence, removable_region_boundaries)
        for genome in genomes:
            genome_blast_result = blast_against_genome(purged_sequence, genome)
            print(json.dumps(genome_blast_result, indent=2))


def blast_against_genome(sequence, genome):
    db_filepath = get_db_filepath(genome, "nhr")
    temp_filepath = write_sequence_to_temp_file(sequence)
    blast_command = NcbiblastnCommandline(cmd="tblastn",
                                          outfmt=15,
                                          query=temp_filepath,
                                          db=db_filepath)
    return json.loads(blast_command()[0])


def replace_string_segment(start, end, string, replacement_char="X"):
    beginning = string[:start]
    ending = string[end+1:]
    distance = (end - start) + 1
    return beginning + replacement_char * distance + ending


def remove_regions_from_sequence(sequence, regions):
    sequence_string = sequence.sequence
    for start, end in regions:
        sequence_string = replace_string_segment(start, end, sequence_string)
    sequence.sequence = sequence_string
    return sequence


def boundaries_match(hit_one, hit_two):
    return (hit_one["query_from"] == hit_two["hit_from"] and hit_one["query_to"] == hit_two["hit_to"]) or (
                hit_two["query_from"] == hit_one["hit_from"] and hit_two["query_to"] == hit_one["hit_to"])


def bracketed_region_information(hit):
    query_from = hit["query_from"]
    query_to = hit["query_to"]
    hit_from = hit["hit_from"]
    hit_to = hit["hit_to"]
    if query_from < hit_from:
        return {
            "size": hit_from - query_to,
            "start": query_from,
            "end": hit_to
        }
    else:
        return {
            "size": query_from - hit_to,
            "start": hit_from,
            "end": query_to
        }


def deletable_region(region):
    if "size" in region:
        return region["start"], region["end"]
    elif "evalue" in region:
        return region["query_from"], region["query_to"]


def find_regions_to_remove(blast_data):
    blast_hits = blast_data["BlastOutput2"][0]["report"]["results"]["bl2seq"][0]["hits"][0]["hsps"]
    blast_score_groups = {}
    for hit in blast_hits:
        if (hit["bit_score"], hit["evalue"]) not in blast_score_groups:
            blast_score_groups[(hit["bit_score"], hit["evalue"])] = [hit]
        else:
            blast_score_groups[(hit["bit_score"], hit["evalue"])].append(hit)
    ltr_regions = []
    non_ltr_repeats = []
    for key, hits in blast_score_groups.items():
        if len(hits) != 2:
            continue
        else:
            if boundaries_match(hits[0], hits[1]):
                bracketed_region = bracketed_region_information(hits[0])
                if bracketed_region["size"] < 500:
                    non_ltr_repeats += hits
                else:
                    ltr_regions.append(bracketed_region)
    return ltr_regions + non_ltr_repeats


def write_sequence_to_temp_file(sequence):
    filepath = Path(f"/tmp/{sequence.start}-{sequence.end}_{sequence.accession_id}_{sequence.genome}.fasta")
    with open(filepath, 'w') as temp_outfile:
        temp_outfile.write(f">{filepath.stem}\n")
        temp_outfile.write(sequence.sequence)
    return str(filepath)


def blast_sequence_against_itself(sequence):
    temp_filepath = write_sequence_to_temp_file(sequence)
    blast_command = NcbiblastnCommandline(query=temp_filepath,
                                          subject=temp_filepath,
                                          outfmt=15,
                                          strand="both")
    # Return is a tuple, we only want the actual data returned
    return json.loads(blast_command()[0])


def expand_sequence_ranges(start, end, expansion_size):
    if start - expansion_size >= 0:
        return start - expansion_size, end + expansion_size
    elif start - (expansion_size / 2) >= 0:
        LOGGER.info("Expansion size exceeded length of sequence, expanding instead by half the amount")
        expansion = expansion_size / 2
        return start - expansion, end + expansion
    else:
        LOGGER.warning("Cannot expand sequence due to start point limitations")
        return start, end


def get_db_filepath(genome, filetype):
    db_dir = Path(homify_path(CONFIG.genome_db_storage))
    for file in db_dir.glob(f"*.{filetype}"):
        if file.stem == genome:
            if filetype in ["nhr", "nin"]:
                return str(file.parent) + "/" + file.stem
            else:
                return str(file)
    raise GenomeDbNotFoundException(f"Cannot find fasta file for genome {genome}")


def get_from_genome_db(start, end, sequence):
    genome_db_filepath = get_db_filepath(sequence.genome, "fasta")
    for seq_record in SeqIO.parse(genome_db_filepath, "fasta"):
        if seq_record.id == sequence.accession_id:
            return str(seq_record.seq[start:end])
    raise AccessionIdNotFoundException(f"Cannot find record for accession id {sequence.accession_id}"
                                       f" in genome {sequence.genome}")


def expand_sequence(sequence, expansion_size=5000):
    expanded_start, expanded_end = expand_sequence_ranges(sequence.start, sequence.end, expansion_size)
    expanded_sequence = get_from_genome_db(expanded_start, expanded_end, sequence)
    return Sequence(expanded_start, expanded_end, expanded_sequence,
                    sequence.accession_id, sequence.genome)