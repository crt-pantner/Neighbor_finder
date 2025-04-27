import csv
import logging
from pathlib import Path

from Bio import SeqIO

# Handles opening the fasta file
def open_fasta_file(location):
    """Helper function that simply opens the fasta files and loads the into a seqio dictionary"""

    logging.debug(f"Opening fasta file: {location}")
    dict = SeqIO.to_dict(SeqIO.parse(location, "fasta"))
    return dict

def output_pairs(protein, pairs, output_file):
    output = []

    output_file = Path(output_file)
    output_file = output_file.resolve() / "pairs"

    # Preping query protein for output
    output.extend(
        [protein.short_name, protein.protid, protein.feature.start, protein.feature.end]
    )

    # Preping the pairs for output
    for f in pairs:
        split_id = f.id.split("_")
        short_name = split_id[0]
        protid = split_id[1]
        output.extend([short_name, protid, f.start, f.stop])

    # outputing to file.
    with open(f"{str(output_file)}.csv", "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(output)