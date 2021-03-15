import click
from click import echo
import logging

import blast_probes

DEFAULT_ALIGNMENT_LENGTH_THRESHOLD = 400
DEFAULT_E_VALUE_THRESHOLD = 0.009

LOGGER = logging.getLogger("ervin")
LOGGER.setLevel(logging.INFO)


@click.group()
def main():
    pass


@main.command("blast-probes")
@click.option("-f", "--filename", required=True, type=click.Path(), help="Input fasta file path")
@click.option("-gdb", "--genome-db", required=True, type=str,
              help="Name of the genome database against which to BLAST the input sequences")
@click.option("-o", "--output-dir", type=click.Path(), default=blast_probes.default_output_dir(),
              help="Directory to which the output should be written")
@click.option("-a", "-alignment-threshold", "alignment_threshold", type=int, default=DEFAULT_ALIGNMENT_LENGTH_THRESHOLD,
              help="Minimum length threshold that BLAST result alignment sequence lengths should exceed")
@click.option("-e", "--e-value-limit", type=float, default=DEFAULT_E_VALUE_THRESHOLD,
              help="Maximum e-value by which to threshold the results returned by the BLAST run")
def blast_probes_cli(filename, genome_db, output_dir, alignment_threshold, e_value_limit):
    blast_probes.run_blast_probes(filename, genome_db, output_dir, alignment_threshold, e_value_limit)


if __name__ == "__main__":
    main()
