Ta skripta (test_gffutils.py) ti poišče zaželen gen ven iz gff datoteke. 

Small utility script that allows matching of neighboring proteins.

As an input it accepts a query fasta file from JGI, in which there are proteins whose protein neighbourhoods you want to explore. 
It is also necessary to provide a folder which contains gff files belonging to the proteins.

For each protein it's corresponding gff file is found. Using gffutils a sqlite database for the protein is created in memory, and the feature for the protein gets extracted. Currently, only genes get extracted.

The last necessary input to provide is the fasta file of proteins around query proteins you want to search.

By deafault, the whole scaffold of the query protein gets extracted, and the program checks wether target protein ids are present on the scaffold. If they are the query and it's pair gets output to an outfile.

By using the optional -choords flag a specific flanking region can be set. For example -coords 10 000 will find proteins located 10 000 bp upstram and downstream from the start and end of the query protein.

For the species files, whose gff3 files can not be found, the species get output to a not_found.txt file.

15. 05. 2025: 
    poskusil jo bom uporabiti za namene gledanja kateri geni so sosednji

Treba je izboljšat output za pare. trenutno je messy in brez headerja. In se vedno naloži v isto datoteko.

Treba je dodati funkcionalnost, s katero lahko obdržiš nastale db datoteke, namesto da se hranijo v memory. To ni striktno nujno ampak je uporabno, če delaš primerjave.

Preverjeno na aegerolizinih in 10 000 bp levo in desno najde vedno ustrezen par, ki je vedno macpf. Sedaj sem 99.99% prepričan, da bo pare vedno pravilno našlo.

