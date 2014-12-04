Quick Information
total snp to omim mappings (in snp_omim): 1339
total rs numbers (in 17E.txt): 6521408
total records in (omim.txt 2): 23836

********************************************************************************

.ipynb Files:
These files can be run inside of ipython notebook if you have the anaconda
distribution of python installed. If you want to import a .ipynb file into
another .ipynb file I would suggest downloading it as a .py file first 
instead of the file tab for the ipython notebook application. 
Note - The .ipynb files will be downloaded into .py files after every 
use on the branch promethease so if you don't want to deal with .ipynb 
then you don't have too. 

********************************************************************************

Information File Descriptions:
17E.vcf - 
Contains information on individual X and the variations that afflict this 
individual. 

cytoBand.txt - 
Cytogenetic information along with the interval positions that correspond.
Obtained from hg19.

genemap - 
Information linking omim records with the corresponding cytogenetic location.
Obtained from the omim database.

genemap2.txt 
Like gene map.
Obtained from the omim database.

genemap.key -
Key describing the information inside of the genemap and morbid map file.
Obtained from the omim database.

morbidmap -
Condensed disease information linked to a cytogenetic location. 
Obtained from the omim database.

mim2gene.txt -
Links the mim id to some gene id.  
Obtained from the omim database.

omim.txt -
The main omim database file. Contains descriptions of diseases and many
other fields of information.  

snp_omim -
Used to map from the rs id to the mim id. 
Obtained from the dbSNP database.

********************************************************************************

Python File Descriptions:
omim_parsing.py/omim_parsing.ipynb - 
Extracts information from omim database files. Currently only parses the
omim.txt file for each of the records.

promethease.ipynb - 
Imports omim_parsing.py as well as vcf_parsing.py and uses their functions
to currently map from the vcf file 

vcf_parsing.py/vcf_parsing.ipynb - 
Parses the 17E.vcf into rs id's. 

 
