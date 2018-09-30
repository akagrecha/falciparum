import re , operator, random, pdb
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord

_start = r'ATG'
_stop = r'(?:TAG|TGA|TAA)'
_nonstop = r'(?:[CAG][TCAG]{2}|T(?:[TC][TCAG]|[AG][TC])|TGG)'
_codon = r'(?:[TCAG]{3})'
# re.I to ignore case
_orf_re = re.compile('(ATG'+ _nonstop + '*)\
('+ _stop +')',flags = re.I) # Edited for finding only orfs
_lead_re = re.compile(r'ATG',flags = re.I)


def make_id(in_id):
    id_re = re.compile("(^[\w\d]+\.?\d?)")
    id = id_re.findall(in_id)[0]
    return id

def extract_orfs(record):
    """
    :param record: Bio.SeqRecord.SeqRecord
        record extracted from parsed input
    :return orfs: list
        open reading frame in the genome
    """
    orf_id = make_id(record.id)
    orfs = []
    ids = []
    start_pos = []
    seqt = record.seq
    # finding all starting positions (ATG) in the genome
    for a in _lead_re.finditer(str(seqt)):
        start_pos.append(a.start())
    for i in start_pos:
        a = _orf_re.search(str(seqt),i)
        if a is not None:
            s = seqt[int(a.start()-utr):int(a.end()+dtr)]
            if len(s) >= int(min_len+utr+dtr) and len(s) <=int(max_len+utr+dtr) and len(s[utr:-dtr])%3 == 0:
                new_orf_id = "%s:%s..%s(+) %s:0..%s(1) Upstream of %s"%\
                             (orf_id,int(a.start()),int(a.end()),orf_id,len(seqt),orf_id)
                new_orf = SeqRecord(s,id = new_orf_id, description = record.description)
                if new_orf.id not in ids:
                    orfs.append(new_orf)
                    ids.append(new_orf.id)
    return (orfs)

utr = 0
dtr = 0
min_len = 1
max_len = 100000000

import sys
assert(len(sys.argv) == 3)
file_handle = sys.argv[1]
output_file = sys.argv[2]

newlist = []

for record in SeqIO.parse(file_handle,"fasta"):
    for i in extract_orfs(record):
        newlist.append(i)

SeqIO.write(newlist, output_file, "fasta")
