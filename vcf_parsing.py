
# coding: utf-8

# In[3]:

__author__ = 'joaquinreyna'

import pandas as pd
import math

def find_interval(vcf_pos, interval_list):
    """Given a list of a chromosome find the interval where the
    vcf position is found.
    Argument(s):
        vcf_pos: The position where the current vcf entry is found.
        interval_list: The intervals for the chromosome bands.
    
    Return(s):
        Tuple with the interval in which the vcf position is found.
    """
    
    for lower_bound, upper_bound in interval_list:
        if lower_bound <= vcf_pos and vcf_pos <= upper_bound:
            return (lower_bound, upper_bound)
    return "Roger we have a problem."

def create_rs_omim_df(file_n = "snp_omim"):
    """Given the snp_omim file from the dbSNP website a dataframe 
        is made storing rs and mim ids. 
        Argument(s):
            The current working snp_omim file.
        Return(s): 
            Dataframe - Indices: rs id, Columns: mim id
    """
    
    with open(file_n) as so_mapping:
        so_dict = dict()
        for mapping in so_mapping:
            mapping = mapping.strip().split()
            rs_id, mim_id = mapping[0:2]
            so_dict[rs_id] = mim_id
        index = ["{}".format(i) for i in range(len(so_dict))]
        so_series = pd.Series(so_dict)
        so_df = pd.DataFrame(so_series, columns = ["mim_id"])
        so_df.sort(columns = "mim_id", axis = 0, inplace = True)
    return so_df

class vcf():
    """Reads and stores vcf information."""
    
    def __init__(self, file_n = "17E.vcf"):
        self.file_name = file_n

    def create_cytogenetic_dict(self):
        """Making a dictionary for the cytogenetic locations."""
        
        with open("cytoBand.txt") as cyto_file:
            #initiating the cyto_dictionary
            self.cyto_dict = dict()
            initiate = cyto_file.next().strip().split()
            prev_chrom = initiate[0].replace("chr","")
            for pos in range(1,3):
                initiate[pos] = int(initiate[pos])
            prev_pos = tuple(initiate[1:3])
            self.cyto_dict[prev_chrom] = dict()
            self.cyto_dict[prev_chrom][prev_pos] = initiate[3]

            #Reading the cytoBand.txt file and filling cyto_dict
            for cyto_loc in cyto_file:
                cyto_loc = cyto_loc.strip().split()
                cyto_loc[0] = cyto_loc[0].replace("chr", "")
                
                for pos in range(1,3):
                    cyto_loc[pos] = int(cyto_loc[pos])
                cur_chrom = cyto_loc[0]
                cur_pos = tuple(cyto_loc[1:3])
                
                if cur_chrom == prev_chrom:
                    self.cyto_dict[cur_chrom][cur_pos] = cyto_loc[3]

                else:
                    self.cyto_dict[cur_chrom] = dict()
                    self.cyto_dict[cur_chrom][cur_pos] = cyto_loc[3]
                    prev_chrom = cur_chrom

    def create_vcf_interval_dict(self):
        """Opening the vfc file and reading the records for the genes affected.
            Return(s):
                Dictionary - Key(s): position interval, Value(s): cytogenetic position 
        """
        
        create_cytogenetic_dict()
        with open(self.file_name) as vcf_file:
            vcf_cytohits_dict = dict()
            vcf_file.next()
            initiate = vcf_file.next().strip().split()
            cur_chrom = initiate[0]
            cur_pos = int(initiate[1])
            vcf_cytohits_dict[cur_chrom] = dict()
            cur_interval = find_interval(cur_pos, self.cyto_dict[cur_chrom].keys())
            vcf_cytohits_dict[cur_chrom][cur_interval] = 1

            for vcf_record in vcf_file:
                vcf_record = vcf_record.strip().split()
                cur_chrom = vcf_record[0].strip()
                cur_pos = int(vcf_record[1])
                cur_interval = find_interval(cur_pos, self.cyto_dict[cur_chrom].keys())

                if cur_chrom in vcf_cytohits_dict.keys():
                    if cur_interval in vcf_cytohits_dict[cur_chrom].keys():
                        vcf_cytohits_dict[cur_chrom][cur_interval] += 1
                else:
                    vcf_cytohits_dict[cur_chrom] = dict()
                    vcf_cytohits_dict[cur_chrom][cur_interval] = 1


    def create_vcf_rs_list(self):
        """Reading the vfc file and to obtain a list of all the rs id's."""
        with open(self.file_name) as vcf_file:
            vcf_rs_ids = list()
            vcf_file.next()
            
            for vcf_record in vcf_file:
                vcf_record = vcf_record.strip().split()
                if vcf_record[2].startswith("rs"):
                    vcf_rs_ids.append(vcf_record[2].replace("rs", ""))
                    
        return vcf_rs_ids


# In[ ]:



