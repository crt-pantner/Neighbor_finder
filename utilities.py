from pathlib import Path
import gffutils
import logging

from helpers import transform
from protein_class import Protein

def database_creator(db_path, gff_file_path, force=False):

    
    gff_file_path = Path(gff_file_path)
    
    

    gffutils_database = gffutils.create_db(
        str(gff_file_path.resolve()),
        dbfn=str(db_path),
        merge_strategy="create_unique",
        id_spec={"transcript": "transcriptId", "gene": "ID", "exon": "ID", "mRNA": "ID"},
        gtf_transcript_key="transcriptId",
        gtf_gene_key="gene",
        transform=transform,
        force=force
    )
    return gffutils_database


def get_organism_objs(seqio_dict):
    """Creates instances of organism objects form the seqio dictionary. This just ensures easier handling downstream"""

    organism_objs = []

    for key in seqio_dict.keys():
        org_name = key.split("|")[1]
        protein_id = key.split("|")[2]
        organism = Protein(
            short_name=org_name,
            long_name=key,
            sequence=seqio_dict[key],
            protid=protein_id,
        )
        organism_objs.append(organism)
    return organism_objs

def get_logger(args):
    level = logging.basicConfig(
        level=logging.DEBUG if args == True else logging.WARNING, format="%(message)s"
    )

def find_gff_file(organism, outpath, gff_folder=None):
    """Finds the gff file for each organism."""

    logging.debug(f"\nLooking for {organism} gff files")
    
    paths = [
    p
    for p in Path(gff_folder).rglob(f"{organism}*.gff3*")
    if not (p.name.endswith(".tar") or p.name.endswith(".tar.gz"))
]
    

    if paths:
        logging.debug("Paths of gff3 files found.")
    else:
        logging.debug(f"Gff file for {organism} not found, writing to output.")
        with open(f"{outpath}_not_found.txt", "a", newline="\n") as outfile:
            outfile.writelines(f"No gff file found: {organism}\n")
        raise FileNotFoundError

    return [str(p) for p in paths]