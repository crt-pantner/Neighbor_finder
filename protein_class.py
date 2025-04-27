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