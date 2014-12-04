
# coding: utf-8

# In[12]:

import pandas as pd
import omim_parsing as op
import vcf_parsing as vp


## Merge: VCF to SNP Omim

# In[13]:

vcf_obj = vp.vcf()
rs_list = vcf_obj.create_vcf_rs_list()


# In[14]:

vcf_df = pd.DataFrame()
vcf_df["rs_id"] = rs_list


# In[15]:

snp_omim_df = vp.create_rs_omim_df()
vcf_snp_df = pd.merge(left = vcf_df, right = snp_omim_df,                      left_on = "rs_id", right_index = True, how = "right")
vcf_snp_df.dropna(inplace = True, subset = ["mim_id"])


## Merge: VCF_SNP_OMIM to Omim DB

# In[16]:

omim_obj = op.omim()
omim_dict = omim_obj.create_omim_dict()


# In[17]:

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


# In[18]:

vcf_omim_df = pd.merge(left = vcf_snp_df, right = omim_df, left_on = "mim_id",                      right_on = "mim_id", how = "inner")
vcf_omim_df.to_csv("vcf_omim_df.csv")


## Grouping By MIM ID

# In[19]:

mim_id_grouping = vcf_omim_df.groupby("mim_id")


# In[ ]:

omim_variation_hits = mim_id_grouping["rs_id"].agg("count")
omim_variation_hits.sort(ascending = False)
omim_variation_hits.to_csv("omim_variation_hits.csv")


# In[ ]:

omim_rs_id_hits = mim_id_grouping["rs_id"].agg(lambda x: x)
omim_rs_id_hits.to_csv("omim_rs_id_hits.csv")


# In[ ]:



