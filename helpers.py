

def get_nighbourhood(start, end, choords):
    """get the desired genomic neighbourhood"""

    upstream = start - choords
    downstream = end + choords
    return upstream, downstream

def transform(d):
    """Helper function for gffutils, ensures that features get loaded in properly."""
    try:

        d["CDS"] = d["exon"].replace("exon_", "CDS:")

    except KeyError:

        pass

    return d