"""
OMIM (Online Mendelian Inheritance in Man) parser.
"""
# Python imports
import re

# Ibidas imports
from parsers import liner, utils
from container import table, scalar

def _splitAndClean(s):
    vars = [x.strip() for x in s.split(',')]
    return vars


class OMIMFlusher:

    OMIM_STATUS_NAMES= [
            'status confirmed', 'status provisional', 
            'status inconsistent', 'status limbo'
    ]

    OMIM_STATUS_CODES = ['C', 'P', 'I', 'L']
    
    OMIM_METHOD_CODES = [
            'A', 'AAS', 'C', 'Ch', 'D', 'EM', 'F', 'H', 'HS', 'L', 'LD', 'M', 
            'OT', 'Pcm', 'Psh', 'R', 'RE', 'REa', 'REb', 'REc', 'REf', 'REl',
            'REn', 'S', 'T', 'V', 'X/A'
    ]
    
    def __init__(self, cd, omim, morbidmap=None, genemap=None):
        self.cd = cd
        self.omim = omim
        self.morbidmap = morbidmap
        self.genemap = genemap
        self._getPrerequisites()

    def _getPrerequisites(self):
        """Retrieves the terms needed for flushing the info from the database."""
        cd = self.cd
        method_codes = ['OMIM:' + x for x in self.OMIM_METHOD_CODES]
        status_names = self.OMIM_STATUS_NAMES
        terms = status_names + method_codes
        self.terms = cd.term[_.name.within(*terms), ('term_id', 'name')]
        assert len(self.terms) == len(terms), \
                "Not all OMIM symbols could be found in the database."
        codes = self.OMIM_STATUS_CODES + self.OMIM_METHOD_CODES
        self.terms = self.terms %+ table(codes, ['OMIM_code'])

    def flush(self):
        cd = self.cd
        self.omim_ds_id = utils.resolveSet(cd, 'OMIM')
        self._flushOMIM()
        #self._flushGenemap()
        cd.term.flush()
        cd.commit()

    def _flushOMIM(self):
        cd = self.cd

        # get items from omim
        omim_recs = self.omim.cont
        omim_recs = omim_recs.fname(OMIM_ID='identifier')
        omim_recs = omim_recs.select(_, scalar(self.omim_ds_id) / "source_id")
        omim_rec_ids = cd.term.extend(omim_recs)
        self.omim_recs = omim_rec_ids %+ omim_recs

    def _flushGenemap(self):
        cd = self.cd
        item_term = []
        for i in range(len(self.genemap.rec)):
            names = self.genemap.rec[i, ('gene symbols',)]()
            ids = cd.item_prop_accession[_.value.within(*names), ('item_id',)]()

            # check if the ids that we have found point to the same item
            checked = []
            if type(ids) == int:
                checked.append(ids)
            else:
                for id in ids:
                    if not id in checked:
                        checked.append(id)
            
            if len(checked) == 1:
                item_term.append( (checked[0], self.genemap.rec[i, ('OMIM_ID',)] ) )
                
            elif len(checked) == 0:
                print "No gene found for entry: " + str(self.genemap.rec[i]())
            else:
                pass
                # to many records found
        

class OMIMGenemapParser:
    
    GENEMAP_COLS = [
            'numbering', 'month entered', 'day entered', 'year entered',
            'location', 'gene symbols', 'gene status', 'title', 'empty 1', 
            'OMIM_ID', 'method', 'comments', 'empty 2', 'disorders 1',
            'disorders 2', 'disorders 3', 'mouse correlate', 'references'
    ]
    GENEMAP_MAPPER = {
        'numbering': lambda x: float(x),
        'month entered': lambda x: int(x),
        'day entered': lambda x: int(x),
        'year entered': lambda x: int(x),
        'gene symbols': _splitAndClean,
        'method': _splitAndClean,
        'OMIM_ID': lambda x: 'MIM:' + x.strip(),
    }

    def __init__(self, filename):
        self.filename = filename

    def parse(self):
        records = liner.parse(
                self.filename, 
                fieldnames = self.GENEMAP_COLS,
                delimiter = '|',
                mapper = self.GENEMAP_MAPPER,
                hasheader = False,
        )
        self.rec = records



class OMIMMorbidmapParser:

    MORBIDMAP_COLS = ['disorder', 'symbols', 'OMIM_ID', 'location']
    MORBID_MAPPER = {
        'symbols': lambda x: [y.strip() for y in x.split(',')],
        'OMIM_ID': lambda x: 'MIM:' + x,
    }

    def __init__(self, filename):
        self.filename = filename

    def parse(self):
        records = liner.parse(
                self.filename, 
                fieldnames = self.MORBIDMAP_COLS,
                delimiter = '|',
                mapper = self.MORBID_MAPPER,
                hasheader = False,
        )
        self.rec = records

    def flush(self, cd):
        omim_ds_id = utils.resolveSet(cd, 'OMIM') # FIXME: add description
        self.rec = self.rec.fname(OMIM_ID='identifier', disorder='name')
        self.rec = self.rec.select(_, scalar(omim_ds_id) / 'source_id')
        term_ids = cd.term.extend(self.rec)

        cd.term.flush()
        cd.commit()


class OMIMParser:

    FIELD_RE = re.compile("""
        ^\*FIELD\*
        \s+
        (?P<type>\w{2}).*$
        """, flags=re.VERBOSE)

    TI_RE = re.compile("""
        [\*|\^|\#|\%]{0,1}    # prefix
        (?P<nr>)\d+     # number of record?
    """, flags=re.VERBOSE)

    MULTILINE_FIELDS = ['AV', 'TX', 'RF', 'CS', 'ED', 'CN', 'SA']
    OTHER_FIELDS = ['CD']

    def __init__(self, filename):
        self.filename = filename

    def _newRecord(self):
        record = {}
        for field in self.MULTILINE_FIELDS:
            record[field] = []
        record['TX'] = ''
        return record
    
    def _findRecords(self):
        """Reads the OMIM file and devides it into records and fields."""
        handle = open(self.filename)
    
        records = []
        record = self._newRecord()
        cur_record_type = ""
        for line in handle:
            line = line.strip()
            match = self.FIELD_RE.match(line)
            if(line == '*RECORD*'):
                if record.has_key('NO'):
                    if record['TI']:
                        if(not (('REMOVED' in record['TI']) | ('MOVED' in record['TI'][0]) )):
                            records.append(record)
                    cur_record_type = ""
                    record = self._newRecord()
            elif match:
                cur_record_type = match.group('type')
            else:
                if cur_record_type == 'NO':
                    record['NO'] = line.strip()
                elif cur_record_type == 'TI':
                    vars = line.split()
                    match = self.TI_RE.match(vars[0])
                    if match:
                        line = " ".join(vars[1:len(vars)])
                    record['TI'] = line.strip(';').strip()
                elif cur_record_type == 'TX':
                    #record['TX'] = record['TX'] + " " + line.strip()
                    pass
                elif cur_record_type in self.MULTILINE_FIELDS:
                    record[cur_record_type].append(line)
                elif cur_record_type in self.OTHER_FIELDS:
                    record[cur_record_type] = line
                else:
                    print "Found unknown line type (%s)." % cur_record_type

        handle.close()
        self.records = records

    def _toContainer(self):
        ar = []
        for rec in self.records:
            ar.append(
                (
                    'MIM:' + rec['NO'],
                    rec['TI'],
                    rec['TX'].strip(),
                )
            )
        self.cont = table(ar, ['OMIM_ID', 'name', 'description'])

    def parse(self):
        self._findRecords()
        self._toContainer()
