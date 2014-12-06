
# coding: utf-8

# In[1]:

import pandas as pd
import math


## Parsing the VCF File

# In[ ]:

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


def create_vcf_interval_dict(file_n = "17E.vcf"):
    """Opening the vfc file and reading the records for the genes affected.
        Return(s):
            Dictionary - Key(s): position interval, Value(s): cytogenetic position 
    """

    create_cytogenetic_dict()
    with open(file_n) as vcf_file:
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
    with open(file_n = "17E.vcf") as vcf_file:
        vcf_rs_ids = list()
        vcf_file.next()

        for vcf_record in vcf_file:
            vcf_record = vcf_record.strip().split()
            if vcf_record[2].startswith("rs"):
                vcf_rs_ids.append(vcf_record[2].replace("rs", ""))

    return vcf_rs_ids


## Parsing the cytoband.txt File

# In[ ]:

def create_cytogenetic_dict(file_n = "cytoBand.txt"):
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


## Parsing the omim.txt File

# In[2]:

def create_omim_dict(file_n = "omim.txt"):
    """Opening the OMIM.txt.Z file and creating a dictionary.
        Key: omim number, Value: dictionary of key::field, value:: field information."""

    field_types  = ["NO", "TI", "TX", "RF", "CS", 
                    "CN", "CD", "ED", "SA"]
    collected_fields = ["NO", "TI", "TX"]

    with open(file_n) as omim_db:
        omim_dict = dict()
        record_dict = dict()
        field_count = 0
        count = 0
        cur_field = "NO"
        prev_field = "TI"
        prev_line = ""
        omim_db.next()
        omim_db.next()
        cur_record_no = omim_db.next().strip()
        prev_record_no = cur_record_no
        cur_string = ""

        for line in omim_db:

            line = line.strip()

            if "*RECORD*" not in line:

                if any("*FIELD* " + _type in line for _type in field_types):
                    cur_field = line.strip().replace("*FIELD* ", "")
                    if cur_field in collected_fields:
                        record_dict[prev_field] = cur_string
                        cur_string = ""

                        if cur_field == "NO":
                            prev_field = "NO"
                            prev_line = line.strip()

                        elif cur_field == "TI":
                            prev_field = "TI"

                        elif cur_field == "TX":
                            prev_field = "TX"

                    else:
                        cur_field = "NOT COLLECTING"

                elif cur_field == prev_field:

                    if cur_field == "NOT COLLECTING":
                        continue

                    if cur_field == "NO":
                        cur_record_no = line.strip()

                    cur_string = cur_string + " " + line


            elif "*RECORD*" in line:
                omim_dict[cur_record_no] = record_dict
                record_dict = dict()

            prev_line = line

            field_count += 1
            count += 1

    return omim_dict

def create_omim_df(omim_dict):
    NO_list = list()
    TI_list = list()
    TX_list = list()
    #Filling lists from the omim dictionary to
    #make a dataframe. 
    for mim_id, fields in omim_dict.items():
        NO_list.append(mim_id)
        TI_list.append(omim_dict[mim_id]["TI"])
        if "TX" in omim_dict[mim_id].keys():
            TX_list.append(omim_dict[mim_id]["TX"])
        else:
            TX_list.append("NA")

    omim_df = pd.concat([pd.Series(NO_list), pd.Series(TI_list), pd.Series(TX_list)], axis = 1)
    omim_df.columns = ["mim_id", "TI", "TX"]
    return omim_df

def create_description_dict(text):
    """Fills the information of the omim_txt object. 
        Argument(s): 
            text - A string resulting from reading in the omim.txt file.
    """
    descriptions = {}
    description = text.split()
    word_ls = list()
    category = ""
    prev_category = ""
    description_fields = ["DESCRIPTION", "CLINICAL FEATURES", "INHERITANCE", "PATHOGENESIS",                          "CLONING", "GENE FUNCTION", "MAPPING", "MOLECULAR GENETICS",                          "GENETICS", "POPULATION", "NOMENCLATURE", "ANIMAL MODEL", "MODEL"]
    collecting_fields = ["DESCRIPTION", "CLINICAL FEATURES", "GENE FUNCTION"]

    for i, word in enumerate(description):

        if i + 1 < len(description):
            word_plus = word + " " + description[i + 1]

        if word in description_fields or word_plus in description_fields:

            if word in collecting_fields:
                category = word 

            elif word_plus in collecting_fields:
                category = word_plus

            if category != prev_category and i != 0:
                if prev_category in collecting_fields:
                    descriptions[prev_category] = " ".join(word_ls)
                prev_category = category
                word_ls = []

        else:
            word_ls += [word]

    if len(word_ls) != 0 and prev_category in collecting_fields:
        if prev_category in ["CLINICAL FEATURES", "GENE FUNCTION"]:
            word_ls = word_ls[1:]
        descriptions[prev_category] = " ".join(word_ls)
    return descriptions

