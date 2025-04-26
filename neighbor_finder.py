import gffutils
from Bio import SeqIO
import argparse
from pathlib import Path
import csv
from tqdm import tqdm
import logging
import sqlite3


class Protein:
    def __init__(self, short_name, long_name, sequence, protid):
        self.short_name = short_name
        self.path = None
        self.long_name = long_name
        self.sequence = sequence
        self.protid = protid
        self.feature = None

    def __str__(self):
        return f"Protein object: {self.long_name}"


# Get command line arguments
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-iq",
        required=True,
        type=str,
        help="Path to fasta file with proteins around which to search (MACPF fasta file)",
    )
    parser.add_argument(
        "-it",
        required=False,
        type=str,
        help="Path to fastafile with target proteins (Aegerolysin fasta  file).",
    )
    parser.add_argument(
        "-gff",
        required=False,
        type=str,
        help="Path to folder that contains gff files.",
        default="D:/Documents/Magistrska/poskus_celotnega/Mycocosm_download",
    )
    parser.add_argument(
        "-choords",
        required=False,
        type=int,
        help="Flanking region around proteins. If no value is passed, will pick proteins from whole chromosome",
        default=0,
    )
    parser.add_argument(
        "-o",
        required=False,
        type=str,
        help="Prefix for output csv file containing the proteins and their pairs. Default is 'output_pair'",
        default="output_pairs",
    )
    parser.add_argument(
        "-db_folder",
        required=False,
        type=str,
        help="Prefix/path to folder for saving databases. If provided,, will keep the gffutils database file. Usefull if program is ran multiple times, as it speeds up the process. On default will only save database into memory.",
        default=":memory:"
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    args = parser.parse_args()
    return args


# Handles opening the fasta file
def open_fasta_file(location):
    """Helper function that simply opens the fasta files and loads the into a seqio dictionary"""

    logging.debug(f"Opening fasta file: {location}")
    dict = SeqIO.to_dict(SeqIO.parse(location, "fasta"))
    return dict


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


def transform(d):
    """Helper function for gffutils, ensures that features get loaded in properly."""
    try:

        d["CDS"] = d["exon"].replace("exon_", "CDS:")

    except KeyError:

        pass

    return d


def logger(args):
    level = logging.basicConfig(
        level=logging.DEBUG if args == True else logging.WARNING, format="%(message)s"
    )

   
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


def open_dbs(gff_paths, db_path, organism_group):
    """Opens the gffutils sqlite database for each organism."""
    
    dbs = []

    if db_path == ":memory:":
        logging.debug("Creating database and loading into memory")

        # Creates the database and loads into memory
        for fn in gff_paths:
            try:
                db = database_creator(db_path=db_path, gff_file_path=fn)
                dbs.append(db)
            except Exception as e:
                logging.debug(f"Database creation error, skipping")
                continue

        logging.debug(f"Database created and loaded into memory.")

    # If the user chooses an output directory
    else:
        db_path = Path(db_path) 
        db_path = db_path / f"{organism_group}.db"
        

        for fn in gff_paths:
            
            if db_path.exists():
                # Grab the database and instantiate FeatureDB object if the database already exists
                logging.debug(f"Database already exists, reusing.")
                db = gffutils.FeatureDB(str(db_path))
                dbs.append(db)

            else:
                # If the database doesn't exist yet, create it and save into user specified directory.

                try:
                    db_path = db_path.resolve()
                    logging.debug(f"Creating database and saving into {db_path}")
                    for fn in gff_paths:
                        logging.debug(f"Creating database from {fn}")
                        
                        
                        db = database_creator(db_path, gff_file_path=fn)
                        dbs.append(db)
                    logging.debug(f"Database created and saved to: {db_path}")
                except sqlite3.OperationalError: #Happens because sometimes the program finds two of the same gff files for some organisms, based on the name - Example Irplac118 and Irplac1 can return same file, which leads to the program trying to create the same database twice.
                    continue
        


    return dbs


def group_proteins(organisms):
    """Groups proteins by each organism. This ensures that the database for each organism gets created only one time, as it is the rate limiting step."""
    groups = {}
    for protein in organisms:
        org = protein.short_name
        path = protein.path

        if org not in groups:
            logging.debug(f"New group for {org}")
            groups[org] = {"gff_path": path, "proteins": []}

        groups[org]["proteins"].append(protein)

    return groups


def get_feature(protein, dbs, outpath):
    """Finds each protein in the query file in the respective database created for each organism groups gff file."""

    logging.debug(f"Finding the {protein} feature from the database")
    ids = []

    start = None
    end = None
    seq_id = None

    for db in dbs:
        try:
            for feature in db.features_of_type("gene"):

                if "proteinId" in feature.attributes:
                    # out of all the features in the database, we find the protein we are looking for, by the protid.
                    if feature.attributes["proteinId"][0] == protein.protid:
                        logging.debug(f"Feature found!")

                        start = int(feature.start)
                        end = int(feature.end)
                        seq_id = feature.seqid
                        protein.feature = feature

                        if start and end and seq_id:
                            return protein, db

        except Exception as e:
            logging.debug(f"Error in DB: {e}")

        logging.debug(f"Protein {protein.protid} not found in any DB")
        with open(f"{outpath}_not_found.txt", "a") as outfile:
            outfile.writelines(f"Protein was not found in gff file, perhaps file is structure differently: {protein.long_name}\n")
        raise FileNotFoundError
            


def get_pairs(protein, db, aegerolysin_proteins, choords):
    feature = protein.feature

    # If no choordinates are given as command line arguments, find all features on the same chromosome as query protein.
    if choords == 0:
        region = db.region(seqid=feature.seqid, featuretype="gene")
    else:
        protein_start = protein.feature.start
        protein_end = protein.feature.end

        # Get upstream and downstram region based on command line arguments.
        start, end = get_nighbourhood(protein_start, protein_end, choords)
        region = db.region(
            seqid=feature.seqid, start=start, end=end, featuretype="gene"
        )
    try:
        pairs = []

        # For every upstream and downstream "protein" in the region
        for feature in region:
            gene_id = feature.id
            protein_id = feature.attributes["proteinId"][0]

            # Find the proteins that have the same protein id as in the target fasta file.
            """This essentialy checks for each protein in the vicinitiy of query protein, if it's present in the target fasta file. So that would mean it is a target aegerolysin."""
            for key in aegerolysin_proteins:
                key = key.split("|")
                if protein.short_name == key[1] and protein_id == key[2]:
                    pair = feature.id
                    if pair not in pairs:
                        pairs.append(feature)
    except IndexError:
        print(f"No pairs foud for protein {protein}")
    except TypeError:
        print("IDK IT DOES NOT WORK")
        pass
    return pairs


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


def get_nighbourhood(start, end, choords):
    """get the desired genomic neighbourhood"""

    upstream = start - choords
    downstream = end + choords
    return upstream, downstream


def main():
    # Get input arguments
    arguments = get_arguments()

    logger(args=arguments.verbose)

    # Open macpf and aegerolysin fasta files and store them as a dict
    macpf_dict = open_fasta_file(arguments.iq)
    aegerolysin_seq_dict = open_fasta_file(arguments.it)

    # Pack all the proteins into a list, with long and short form of the name as well as the sequence and protein id.
    organisms_obj_list = get_organism_objs(macpf_dict)

    # Group the proteins by species.
    organism_groups = group_proteins(organisms_obj_list)

    # For each protein group (same species), find its gff file
    for group in tqdm(organism_groups, desc="Processing groups"):
        try:
            paths = find_gff_file(group, gff_folder=arguments.gff, outpath=arguments.o)

            # Open the database for gffutils for this file
            dbs = open_dbs(gff_paths=paths, db_path=arguments.db_folder, organism_group=group)

            for protein in organism_groups[group]["proteins"]:

                # Find the protein in the gff file
                protein, db = get_feature(protein, dbs, outpath=arguments.o)

                # Find the proteins that lie around our query protein
                pairs = get_pairs(
                    db=db,
                    aegerolysin_proteins=aegerolysin_seq_dict,
                    protein=protein,
                    choords=arguments.choords,
                )

                if pairs:
                    logging.debug(f"Pairs found for {protein}, writing to output.")

                    # Output the protein and it's pairs to output file.
                    output_pairs(protein=protein, pairs=pairs, output_file=arguments.o)
        except FileNotFoundError:
            continue


if __name__ == "__main__":
    main()
