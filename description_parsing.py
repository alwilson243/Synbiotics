
# coding: utf-8

# In[1]:

import pandas as pd


# In[44]:

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
 

