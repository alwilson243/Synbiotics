
# coding: utf-8

# In[29]:

class omim():
    def __init__(self, file_n = "omim.txt"):
        self.file_name = file_n

    def create_omim_dict(self):
        """Opening the OMIM.txt.Z file and creating a dictionary.
            Key: omim number, Value: dictionary of key::field, value:: field information."""

        field_types  = ["NO", "TI", "TX", "RF", "CS", 
                        "CN", "CD", "ED", "SA"]
        collected_fields = ["NO", "TI", "TX"]

        with open(self.file_name) as omim_db:
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

