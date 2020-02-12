import json

from .ervin_utils import get_config

from Bio.Blast.Applications import NcbiblastnCommandline

from Bio import SeqIO


LEUCINE_SEQUENCES = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]

LRR_BUFFER = 5000


def conserved_domain_hits(domain_fasta, genome_db):
    config = get_config()
    command = NcbiblastnCommandline(
        cmd="tblastn",
        outfmt=15,
        query=domain_fasta,
        db=f"{config.genome_db_storage}{genome_db}",
    )
    query_result = command()
    parsed_resultsets = json.loads(query_result[0])["BlastOutput2"]
    for parsed_result in parsed_resultsets:
        hits = parsed_result["report"]["results"]["search"]["hits"]
        for hit in hits:
            yield hit


def get_accession_id_from_title(title):
    return title.split()[0]


def find_lrr_start(seq):
    find_results = [seq.upper().find(leucine_seq) for leucine_seq in LEUCINE_SEQUENCES]
    actual_results = [result for result in find_results if result != -1]
    if not actual_results:
        return None
    else:
        return min(actual_results)


def run_conserved_domains(filename, genome_db):
    config = get_config()
    hit_dict = {}
    for hit in conserved_domain_hits(filename, genome_db):
        best_hsp = sorted(hit["hsps"], key=lambda i: i["bit_score"], reverse=True)[0]
        hit_from = best_hsp["hit_from"]
        hit_to = best_hsp["hit_to"]
        strand = "positive" if hit_from < hit_to else "negative"
        accession_id = get_accession_id_from_title(hit["description"][0]["title"])

        if strand == "positive":
            hit_coords = (hit_from, hit_to)
        else:
            hit_coords = (hit_to, hit_from)

        if accession_id not in hit_dict:
            hit_dict[accession_id] = [hit_coords]
        else:
            hit_dict[accession_id].append(hit_coords)

    extended_hit_dict = {}
    for count, seq_record in enumerate(
        SeqIO.parse(
            f"/Users/rob/.ervin/genome_db_store/{genome_db}.fasta", format="fasta"
        )
    ):
        if seq_record.id in hit_dict:
            for coords in hit_dict[seq_record.id]:
                seq_str = str(seq_record.seq)
                modified_coords = (coords[0] - LRR_BUFFER, coords[0])
                if modified_coords[0] >= 0:
                    extended_hit_dict[f">{seq_record.id} {coords[0]}:{coords[1]}"] = {
                        "seq": seq_str[coords[0] : coords[1]],
                        "ext": seq_str[coords[0] - LRR_BUFFER : coords[0]],
                    }

    with open("conserved_domain_results.fasta", "w") as fasta_out:
        for key, value in extended_hit_dict.items():
            lrr_start = find_lrr_start(value["ext"])
            if lrr_start is not None:
                fasta_out.write(f"{key}\n{value['ext'][lrr_start:]}{value['seq']}\n")
