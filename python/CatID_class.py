
class CatID():
    """CatID characteristics taken from the Stacks export_sql file."""
    def __init__(self, ID, consensus, num_snps, snps, num_alleles, alleles, deleveraged):
        self.ID = ID
        self.consensus = consensus
        self.num_snps = num_snps
        self.snps = snps,
        self.num_alleles = num_alleles
        self.alleles = alleles
        self.deleveraged = deleveraged