from functools import lru_cache
from pathlib import Path
from dataclasses import dataclass
import json

DEF_DIR = Path(__file__).parent.parent / "virus_definitions"


@dataclass
class VirusRule:

    proteins: dict
    distance_tolerance: int
    combinations: list

    def min_length(self):
        min_length = min(val for key, val in self.proteins.items())
        return min_length * self.distance_tolerance

    def max_length(self):
        return sum(val * self.distance_tolerance for key, val in self.proteins.items())

    def is_protein_sequence_valid(self, virus_def, protein_one, protein_two):
        if True: #protein_two occurs after protein_one in any of the combinations
            for combination in self.combinations:

            return True
        else:
            return False

    def are_proteins_neighbours(self, virus_def, protein_one, protein_two):
        if self.is_protein_sequence_valid(protein_one.protein, protein_two.protein):
            # If distance between proteins is valid (next-door between tolerance, next door+1 within next door plus tolerance, etc)




@lru_cache(maxsize=None)
def get_virus_definitions():
    defs = []
    for filename in DEF_DIR.glob("*.json"):
        with open(filename) as definition:
            defs.append(VirusRule(**json.load(definition)))
    return defs
