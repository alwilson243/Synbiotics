print "do not run this"

file=open("index.html", "wt")

file.write("<ul>\n")
For each record in dataframe

  file.write("<li>\n")
  file.write("Omim Record " + record.omim_id + "\n")
  # http://www.decalage.info/en/python/print_list assuming its an array of ints
  file.write("Associated rs_ids: " + str(record.rs_ids).strip('[]') + "\n")
  file.write("Title: " + record.TI + "\n")
  file.write("Description: " + record.TX["Description"] + "\n")
  file.write("Clinical Features: " + record.TX["Clinical_Features"] + "\n")
  file.write("Gene Function: " + record.TX["Gene_function"] + "\n")
  file.write("</li>\n")

file.write("</ul>\n")