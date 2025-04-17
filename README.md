Ta skripta (test_gffutils.py) ti poišče zaželen gen ven iz gff datoteke. 

15. 05. 2025: 
    poskusil jo bom uporabiti za namene gledanja kateri geni so sosednji

Treba je dodat:
    - funkcionalnost s pomočjo katere lahko notri vpišeš ali bi rad da se ti izloči cel scaffold ali samo flanking zaporedja, in če samo flanking zaporedja potem command line argument s pomočjo katerega določiš širino v eni smeri.

Timer za izračun koliko časa bo trajala celotna operacija, recimo s pomočjo tdm.
    
    pip install tqdm

    for protein in organism_groups[group]["proteins"]:
    start, end, seqid, db = get_feature(protein, dbs)
    ...