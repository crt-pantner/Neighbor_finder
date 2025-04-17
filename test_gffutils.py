import gffutils, re, time
from Bio import SeqIO
import argparse
from pathlib import Path
import csv

class Protein:
    def __init__(self, short_name, long_name, sequence, protid):
        self.short_name = short_name
        self.path = None
        self.long_name = long_name
        self.sequence = sequence
        self.protid = protid
    
    def __str__(self):
        return f"Protein object: {self.long_name}"
        

 #Get command line arguments
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-iq", required=True, type=str, help="Path to fasta file with proteins around which to search (MACPF fasta file)")
    parser.add_argument("-it", required=False, type=str, help="Path to fastafile with target proteins (Aegerolysin fasta  file).")
    parser.add_argument("-gff", required=False, type=str, help="Path to folder that contains gff files.", default="D:/Documents/Magistrska/poskus_celotnega/Mycocosm_download")
    args = parser.parse_args()
    return args

#Handles opening the fasta file
def open_fasta_file(location):
    print(f"Opening fasta file: {location}")
    dict = SeqIO.to_dict(SeqIO.parse(location, "fasta"))
    return dict

def get_organism_objs(seqio_dict):

    organism_objs = []
    
    for key in seqio_dict.keys():
        org_name = key.split("|")[1]
        protein_id = key.split("|")[2]
        organism = Protein(short_name=org_name, long_name=key, sequence=seqio_dict[key], protid=protein_id)
        organism_objs.append(organism)
    return organism_objs

def find_gff_file(organism, gff_folder=None):
    
    print(f"Looking for {organism} gff file")
    paths = list(Path(gff_folder).rglob(f"{organism}*.gff3*"))
    
    if paths:
        print("paths found")
    else:
        print(f"Gff file for {organism} not found, writing to output.")
        with open("no_gff_files.txt") as outfile:
            outfile.writelines(f"{organism}")
        raise FileNotFoundError


        
    return [str(p) for p in paths]
    
def transform(d):

    try:

        d["CDS"] = d['exon'].replace("exon_", "CDS:")

    except KeyError:

        pass

    return d

def open_dbs(paths, group):
    begining_time = time.time()

    print("Creating database and loading into memory")
    
    dbs = []

    for fn in paths:
        try:
            db = gffutils.create_db(fn, dbfn=":memory:", merge_strategy="create_unique",
                                    id_spec={'transcript': 'transcriptId', 'gene': 'Name', "exon":"ID"}, 
                                    gtf_transcript_key='transcriptId', 
                                    gtf_gene_key='gene',
                                    transform=transform)
            dbs.append(db)
        except Exception as e:
            print(f"Database creation error, skipping")
            continue
        
    

    end_time = time.time()

    print(f"Database created and loaded into memory. Elapsed time: {round(end_time - begining_time)} seconds")
    

    organism_names = []
    for organism in organism_names:
        if organism.short_name not in organism_names:
            organism_names.append(organism.short_name)
    else:
        pass

    return dbs


def group_proteins(organisms):
    groups = {}
    for protein in organisms:
        org = protein.short_name
        path = protein.path

        
        
        if org not in groups:
            print(f"New group for {org}")
            groups[org] = {
                "gff_path": path,
                "proteins": []
            }

        
        groups[org]["proteins"].append(protein)

    return groups

def get_feature(protein, dbs):
    print(f"Finding the {protein} feature from the database")
    ids = []

    start = None
    end = None
    seq_id = None

    
    for db in dbs:
        print(db)
        try:
            for feature in db.features_of_type("gene"):
                
                    if "proteinId" in feature.attributes:
                        #out of all the features in the database, we find the protein we are looking for, by the protid.
                        if feature.attributes["proteinId"][0] == protein.protid:
                            print(f"Feature found!")
                            print({feature})
                            start = int(feature.start)
                            end = int(feature.end)
                            seq_id = feature.seqid
                            if start and end and seq_id:
                                return start, end, seq_id, db
            
        except Exception as e:
            print(f"Error in DB: {e}")

        print(f"Protien {protein.protid} not found in any DB")


    

def get_pairs(start, end, seq_id, db, aegerolysin_proteins, protein):
    try:
        pairs = []
        for feature in db.region(seqid=seq_id, start=start-10000, end=end+10000, featuretype="gene"):
            gene_id = (feature.id)
            protein_id = (feature.attributes["proteinId"][0])

            

            #koda ki ti izmed vseh proteinov lahko poišče tistega, ki ima ujemajoč protein id.
            for key in aegerolysin_proteins:
                if protein_id in key:
                    
                    pairs.append(feature.id)
    except IndexError:
        print(f"No pairs foud for protein {protein}")
    except TypeError:
            print("IDK IT DOES NOT WORK")
            pass           
    return pairs
    

def output_pairs(protein, pairs):
    with open("output_pairs.csv", "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([protein.long_name, pairs])

def main():
    #Get input arguments
    arguments = get_arguments()
    
    
    #Open macpf and aegerolysin fasta files and store them into variable
    macpf_dict = open_fasta_file(arguments.iq)

    
    aegerolysin_seq_dict = open_fasta_file(arguments.it)

    #Pack all the proteins into a list, with long and short form of the name as well as the sequence.
    organisms_obj_list = get_organism_objs(macpf_dict)

    #Poišči gff3 datoteko za vsak organizem, če je ni jo vrži ven v neko datoteko
    #organisms_obj_list = find_gff_file(organisms=organisms, gff_folder=arguments.gff)
    
    #grupiraj organizme 
    organism_groups = group_proteins(organisms_obj_list)
    
    for group in organism_groups:
        try:
            paths = find_gff_file(group, gff_folder=arguments.gff)



            dbs = open_dbs(paths, group)

            for protein in organism_groups[group]["proteins"]:
                start, end, seqid, db = get_feature(protein, dbs)
                pairs = get_pairs(start, end, seqid, db, aegerolysin_proteins=aegerolysin_seq_dict, protein=protein)
                if pairs:
                    print(f"Pairs found for {protein}, writing to output.")
                    pairs = ", ".join(pairs)
                    output_pairs(protein=protein, pairs=pairs)
        except FileNotFoundError:
            pass
        


    #organism_groups[group]["proteins"]


        
        


        



    #odpri podatkovno bazo za vsak organizem


    #za vsak protein id poišči okolico proteina

    #za vsako okolico proteina 

main()

"""


    
"""
