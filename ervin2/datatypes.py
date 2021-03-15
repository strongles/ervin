from dataclasses import dataclass
from dataclasses import field


@dataclass
class BlastHit:
    description: list[dict]
    hsps: list[dict]
    num: int
    len: int

    @property
    def desc_dict(self):
        """Get the description for this hit"""
        return self.description[0]

    @property
    def scaffold_id(self):
        """Get the Accession ID for this hit"""
        return self.desc_dict["title"].split()[0]

    @property
    def scaffold_length(self):
        """Get the length of this hit sequence"""
        return self.len

    @property
    def e_value(self):
        """Get the e-value associated with this hit"""
        return self.hsps[0]["evalue"]

    @property
    def alignment_length(self):
        """Get the length of this hit's alignment"""
        return self.hsps[0]["align_len"]

    @property
    def query_sequence(self):
        """Get the query sequence"""
        return self.hsps[0]["qseq"]

    @property
    def hit_sequence(self):
        """Get this hit's sequence"""
        return self.hsps[0]["hseq"]

    @property
    def frame(self):
        """Get the frame of this sequence"""
        return self.hsps[0]["hit_frame"]

    @property
    def start(self):
        """Get the start location of this hit sequence"""
        return min(self.hsps[0]["hit_from"], self.hsps[0]["hit_to"])

    @property
    def end(self):
        """Get the end location of this hit sequence"""
        return max(self.hsps[0]["hit_from"], self.hsps[0]["hit_to"])

    @property
    def direction(self):
        """Get the direction of this hit's sequence"""
        direction_map = {
            True: "P",
            False: "N"
        }
        return direction_map[self.hsps[0]["hit_from"] < self.hsps[0]["hit_to"]]

    def to_tsv_record(self, accession_id):
        tab = '\t'
        if self.direction == "P":
            start = self.start
            end = self.end
        else:
            start = self.end
            end = self.start
        stringified_output = [str(output) for output in
                              [accession_id, self.scaffold_id, self.scaffold_length, start, end, self.e_value,
                               self.alignment_length, self.query_sequence, self.hit_sequence, self.frame]]
        return f"{tab.join(stringified_output)}\n"


@dataclass
class BlastResult:
    query_masking: list[dict]
    hits: list[dict]
    stat: dict
    query_id: str
    query_title: str
    query_len: int

    @property
    def query_accession_id(self):
        """Get the Accession ID for the query sequence"""
        return self.query_title.split()[0]

    @property
    def query_length(self):
        """Get the length of the query sequence"""
        return self.query_len

    @property
    def hit_results(self):
        """Get an iterator over the BLAST hits"""
        for hit in self.hits:
            if len(hit["hsps"]) == 1:
                yield BlastHit(**hit)
            else:
                for i, hsp in enumerate(hit["hsps"]):
                    sub_hit = hit.copy()
                    sub_hit["hsps"] = [hit["hsps"][i]]
                    yield BlastHit(**sub_hit)


@dataclass
class Fasta:
    title: str
    seq: list[str]

    def __str__(self):
        newline = "\n"
        return f"{self.title}\n{newline.join(self.seq)}\n"

    @property
    def sequence(self):
        return "".join(self.seq)
