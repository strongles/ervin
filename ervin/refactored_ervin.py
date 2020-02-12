import click
import sys

from .ervin_utils import DEFAULT_OUTPUT_DIR
from .ervin_utils import format_timestamp_for_filename
from .ervin_utils import ensure_output_dir_exists

from .probe_blaster import run_probe_blaster

from .probe_finder import run_probe_finder

from .virus_blaster import VIRUS_DB_DEFAULT
from .virus_blaster import run_virus_blaster

from .conserved_domains import run_conserved_domains


@click.group(name="ervin")
@click.pass_context
def cli(context):
    pass


def write_summary_file(probe_count, fasta_count, virus_to_counts, output_dir, run_ts):
    destination_dir = ensure_output_dir_exists(output_dir)
    summary_path = destination_dir / f"summary_{run_ts}.txt"
    with open(summary_path, "w") as summary:
        summary.write(f"PROBE HITS FOUND: {probe_count}\n")
        summary.write(f"PROBE HITS POST FILTERING: {fasta_count}\n")
        summary.write("VIRUS HITS BY NAME:\n")
        for virus_name, count in virus_to_counts.items():
            tab_char = "\t"
            summary.write(f"{tab_char}{virus_name}: {count}\n")


@cli.command(name="erv")
@click.option(
    "-f", "--file", help="Input probe file", type=click.File("r"), required=True
)
@click.option(
    "-o",
    "--output_dir",
    help="Location in which to write the output files",
    type=str,
    default=DEFAULT_OUTPUT_DIR,
)
@click.option(
    "-gdb", "--genome_db", help="Reference genome database", type=str, required=True,
)
@click.option(
    "-vdb",
    "--virus_db",
    help="Locally held Virus reference database",
    type=str,
    default=VIRUS_DB_DEFAULT,
)
@click.option(
    "-a",
    "--alignment_len_threshold",
    help="Minimum length threshold that BLAST result alignment sequences should exceed",
    type=int,
    required=False,
)
@click.option(
    "-e",
    "--e_value",
    help="Maximum e-value by which to threshold the results returned by the BLAST run",
    type=float,
    required=False,
)
@click.pass_context
def find_ervs(
    _ctx, file, output_dir, genome_db, virus_db, alignment_len_threshold, e_value
):
    run_ts = format_timestamp_for_filename()
    blasted_probes, probe_count = run_probe_blaster(
        file.name, genome_db, alignment_len_threshold, e_value, run_ts
    )
    probe_finder_fasta, fasta_count = run_probe_finder(
        blasted_probes, alignment_len_threshold, run_ts
    )
    _, virus_to_counts = run_virus_blaster(
        probe_finder_fasta, virus_db, output_dir, run_ts
    )
    write_summary_file(probe_count, fasta_count, virus_to_counts, output_dir, run_ts)


@cli.command(name="domains")
@click.option(
    "-f",
    "--filename",
    help="Fasta file containing the conserved domain to locate",
    type=click.File("r"),
    required=True,
)
@click.option(
    "-gdb",
    "--genome_db",
    help="Genome DB in which to search for the provided conserved domain",
    type=str,
    required=True,
)
@click.pass_context
def conserved_domains(_ctx, filename, genome_db):
    run_conserved_domains(filename.name, genome_db)


if __name__ == "__main__":
    sys.exit(cli())
