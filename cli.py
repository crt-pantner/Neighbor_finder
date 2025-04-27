import argparse


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