from sys import argv
from Bio import SeqIO

FILEPATH = "../data/GCA_004026985.1_MyoMyo_v1_BIUU_genomic.fna"


def compute_ranges(range_start, range_end):
    return range_start - 5000, range_end + 5000


def get_from_db(range_start, range_end, accession_id):
    # TODO: Find out DB interface
    for seq_record in SeqIO.parse(FILEPATH, "fasta"):
        if seq_record.id == accession_id:
            return seq_record.seq[range_start:range_end]
    raise Exception("Accession ID not found.")


if __name__ == "__main__":
    acc_id = argv[1]
    start = int(argv[2])
    end = int(argv[3])
    # acc_id = "PVIZ010000006.1"

    new_range_start, new_range_end = compute_ranges(start, end)
    result = get_from_db(new_range_start, new_range_end, acc_id)

