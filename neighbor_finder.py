import gffutils
from Bio import SeqIO
import argparse
from pathlib import Path
import csv
from tqdm import tqdm
import logging
import sqlite3

from cli import get_arguments
from utilities import *
from protein_class import Protein
from input_otuput import *
from helpers import get_nighbourhood








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


def get_feature(protein: Protein, dbs, outpath):
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








def main():
    # Get input arguments
    arguments = get_arguments()

    get_logger(args=arguments.verbose)

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
