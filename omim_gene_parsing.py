__author__ = 'joaquinreyna'

class omim_genes():
    def __init__(self, file_n = "mim2gene.txt"):
        self.filename = file_n
    def parse(self):
        """Opening the gene map (from the omim database) and creating a dictionary.
            Key: cytogenetic location and Value: associated mim number(s)."""
        with open("genemap") as omim_gene_map:
            omim_dict = dict()
            for record in omim_gene_map:
                record = record.split("|")
                if record[4] not in omim_dict:
                    omim_dict[record[4]] = [record[9]]
                else:
                    omim_dict[record[4]].append(record[9])

            return omim_dict