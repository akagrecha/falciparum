import numpy as np
from Bio import SeqIO
import csv, re, pickle
from itertools import product

###### Reading and House Keeping ##########
_cfDict = pickle.load(open('cfDict','r'))

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def _generate_inputs(orf_fasta, genome_fasta):
    """
    FUNCTION SIGNATURE
	------------------
	param orf_fasta: str
        path to file containing the file containing all the orfs
    param genome_fasta: str
        path to file containing the file containing all the constructs
    return orf_pos: dict
        key:- ORF, value:- (start_pos, end_pos)
    return orf_dict_arranged:
        key:- construct , value: - ORFs in the construct sorted according to their start position
    return orf_fourth: dict
        key:- ORF, value:- base at 4th position
    return pcomp_dict: dict
        key:- ORF, value:- pcomp of ORF
    return genome_ids: list
        list of constructs
    return input_cdss: list
        list of input coding sequences
    """
    # orf_dict: dict with keys as construct and values as ORFs
    orf_dict = _get_orfs_in_genome(orf_fasta,genome_fasta)
    # orf_pos:
    orf_pos = _get_orf_positions(orf_fasta)

    orf_dict_arranged = _get_arranged_dict(orf_pos, orf_dict)
    construct_dict = {}
    for key in orf_dict:
        orf_list = orf_dict[key]
        for orf in orf_list:
            construct_dict[orf] = key

    construct_max_len = dict(zip(orf_dict.keys(), [0 for i in range(len(orf_dict.keys()))]))
    construct_max_id = dict(zip(orf_dict.keys(), ['' for i in range(len(orf_dict.keys()))]))
    orf_fourth = {} # with keys as ORFs and value as +4 position
    pcomp_dict = {}
    for record in SeqIO.parse(orf_fasta,"fasta"):
        orf = record.id
        construct = construct_dict[orf]
        orf_fourth[orf] = record.seq[4]
        string = str(record.seq)
        l = len(string)
        if l > construct_max_len[construct]:
            construct_max_len[construct] = l
            construct_max_id[construct] = orf
        pcomp = _codon_adap_index(string)
        pcomp_dict[record.id] = pcomp

    input_cdss = []
    input_cdss = construct_max_id.values()
    genome_ids = orf_dict.keys()

    return  orf_pos, orf_dict_arranged, orf_fourth, pcomp_dict, genome_ids, input_cdss

def _get_seq(fasta):
    return [i for i in SeqIO.parse(fasta,"fasta")]

def _get_orfs_in_genome(orf_fasta, genome_fasta):
    """
    FUNCTION SIGNATURE
	------------------
	param orf_fasta: str
        path to file containing the file containing all the orfs
    param genome_fasta: str
        path to file containing the file containing all the constructs
    return upstream_main_gene_dict: dict
        key:-construct, value:- list of orfs in the genome with the construct
    """
    upstream_regec = re.compile("Upstream\sof\s([\w\d]*\.?\d?)")
    upstream_main_gene_dict = {}
    for i in _get_seq(genome_fasta):
        upstream_main_gene_dict.update({i.id:[]})
    for i in _get_seq(orf_fasta):
        try:
            main_gene = upstream_regec.findall(i.description)[0]
            if i.id not in upstream_main_gene_dict[main_gene]:
                upstream_main_gene_dict[main_gene].append(i.id)
        except:
            pass
    return upstream_main_gene_dict

def _get_orf_positions(orf_fasta):
    """
    FUNCTION SIGNATURE
	------------------
    param orf_fasta: str
        path to file containing the file containing all the orfs
    return positions: dict
        key:- ORF, value:- (start_position of ORF, end_position of ORF)
    """
    regec = re.compile("\w*\:(\d*)\.\.(\d*)\(([\w\d\+\-]{1,})")
    positions = {}
    for i in SeqIO.parse(orf_fasta,"fasta"):
        desc = regec.findall(i.description)
        if desc[1][2] == "1":
            start_pos = int(desc[0][0]) - int(desc[1][0])
            end_pos = int(desc[0][1]) - int(desc[1][0])
            positions.update({i.id:(start_pos,end_pos)})#Checked#strt pos
        elif desc[1][2] == "-1":
            start_pos = int(desc[1][1]) - int(desc[0][0])
            end_pos = int(desc[1][1]) - int(desc[0][1])
            positions.update({i.id:(start_pos,end_pos)})# double checked#strt pos
    return positions

def _arrange_according_to_pos(orf_list, orf_pos):
    """
    FUNCTION SIGNATURE
	------------------
    param orf_file: list
        list containing the relevant ORFs
    param orf_pos: dict
        key:- ORFs , value:- (start_position of ORF, end_position of ORF)
    return l: list
        list of ORFs sorted by their starting position
    """
    orf_with_pos = [[i,orf_pos[i][0]] for i in orf_list]
    orf_with_pos.sort(key=lambda x: x[1])
    l = [i[0] for i in orf_with_pos]
    return l

def _get_arranged_dict(orf_pos, orf_dict):
    """
    FUNCTION SIGNATURE
	------------------
    param orf_pos: dict
        key:- ORF, value:- (start_position of ORF, end_position of ORF)
    param orf_dict: dict
        key:- construct, value:- ORFs in the construct
    return orf_dict_arranged: dict
        key:- construct , value: - ORFs in the construct sorted according to their start position
    """
    orf_dict_arranged = {}
    for i,j in orf_dict.items():
        orf_dict_arranged.update({i:_arrange_according_to_pos(j,orf_pos)})
    return orf_dict_arranged

############ preach calculation ##############

def _codon_adap_index(string):
    """
    FUNCTION SIGNATURE
    ------------------
    param string: str
        sequence of ORF
    return cai: float
        codon adaptation index
    """
    if(len(string) == 6):
        return 0
    assert(len(string)%3 == 0)
    assert(len(string)>6)
    assert(string[:3] == 'ATG')
    stopCodons = ['TAA','TAG','TGA']
    assert(string[-3:] in stopCodons)
    codons = [string[i:i+3] for i in range(3,len(string)-3,3)]
    weights = np.array(map(_cfDict.get, codons))
    numCodons = (len(string)-6)/3
    cai = np.prod(weights)**(1.0/numCodons)
    return cai

def _reinit_func(n):
    """
    FUNCTION SIGNATURE
    ------------------
    param n: int
        intercistronic length
    return _: float
        probability of reinitiating
    """
    # preinit is nearly 0.25 at length 25
    if n < 0:
        return 0
    return (np.exp(n/50)-1)/(np.exp(n/50)+1)

def _nd_func(n):
    """
    FUNCTION SIGNATURE
    ------------------
    param n: int
        length of orf
    return _: float
        probability of not dissociating
    """
    # at length 90, pnd = 0.95
    # at length 120, pnd = 0.5
    assert(n>=6)
    if(n<=300):
        return 1.0/(1 + np.exp((n-120)/10))
    else:
        return 1.52299e-08 # avoid overflow errors

def _ppos_func(val):
    """
    FUNCTION SIGNATURE
    ------------------
    param val: str
        base at +4 position
    return _: float
        probability of initiating
    """
    if val == 'A':
        return 0.139
    elif val == 'T':
        return 0.562
    elif val == 'G':
        return 0.729
    elif val == 'C':
        return 0.364
    else: return

def find_total_preinit(n, starts, ends, path):
    """
    FUNCTION SIGNATURE
    ------------------
    param n: int
        index of the CDS
    param starts: np.array(dtype=int)
        start positions of ORFs
    param ends: np.array(dtype=int)
        end positions of ORFs
    param path: np.array(dtype=float)
        list of ones and zeros, one -> initiated in the path
    return p: float
        product of different preinits
    """
    assert(len(starts) == n+1)
    assert(len(ends) == n+1)
    oneInds = np.where(path==1)[0]
    oneInds = np.append(oneInds, n)
    p = 1
    if len(oneInds) == 1:
        return p
    else:
        for i in range(len(oneInds)-2, -1, -1): # after the first one, before the CDS
            el = oneInds[i]
            eln = oneInds[i+1]
            p *= np.multiply.reduce([1 - _reinit_func(starts[j]-ends[el]) for j in range(el, eln)])
            p *= _reinit_func(starts[eln] - ends[el])
        return p

# find the number of uORFs
# find pinits, pelons, pnds, preinits
# def find_ptran(n, pinits, pelons, pnds, starts, ends):
#     """
#     FUNCTION SIGNATURE
#     ------------------
#     param n: int
#         index of CDS
#     param pinits: np.array(dtype=float)
#         pinits of ORFs
#     param pelons: np.array(dtype=float)
#         pelons of ORFs
#     param pnds: np.array(dtype=float)
#         pnds of ORFs
#     param starts: np.array(dtype=int)
#         start positions of ORFs
#     param ends: np.array(dtype=int)
#         end positions of ORFs
#     return ptran: float
#         probability of translation of CDS
#     return preach: float
#         proabability of reaching the CDS
#     """
#     pinits = np.array(pinits)
#     pelons = np.array(pelons)
#     pnds = np.array(pnds)
#     assert(len(pinits) == n+1)
#     assert(len(pelons) == n+1)
#     assert(len(starts) == n+1)
#     paths = np.array(list(product('01',repeat=n)), dtype=float)
#     ptemp_inits = pinits[:n]*pelons[:n]*pnds[:n]
#     ptemp_n_inits = 1 - pinits[:n]
#     ptemps = np.append(ptemp_n_inits.reshape(n,1), ptemp_inits.reshape(n,1))
#     ppaths_temp = [np.product([ptemps[int(el)] for el in pth]) for pth in paths] # erroneous
#     # if a cds is not reached, then 1-pinit should not be multiplied
#     total_preinits = np.array([find_total_preinit(n, starts, ends, path) for path in paths])
#     ppaths = [ppaths_temp[i]*total_preinits[i] for i in range(len(paths))]
#     preach = sum(ppaths)
#     ptran = preach*pinits[n]*pelons[n] # preach*pinit(cds)*pelon(cds)
#     return ptran, preach

def find_ptran(n, pinits, pelons, pnds, starts, ends):
    """
    FUNCTION SIGNATURE
    ------------------
    param n: int
        index of CDS
    param pinits: np.array(dtype=float)
        pinits of ORFs
    param pelons: np.array(dtype=float)
        pelons of ORFs
    param pnds: np.array(dtype=float)
        pnds of ORFs
    param starts: np.array(dtype=int)
        start positions of ORFs
    param ends: np.array(dtype=int)
        end positions of ORFs
    return ptran: float
        probability of translation of CDS
    return preach: float
        proabability of reaching the CDS
    """

    pinits = np.array(pinits)
    pelons = np.array(pelons)
    pnds = np.array(pnds)
    assert(len(pinits) == n+1)
    assert(len(pelons) == n+1)
    assert(len(starts) == n+1)

    pr = np.zeros(n+1)
    pr[0] = 1 # preach for first orf is 1

    ptemp_inits = pinits[:n]*pelons[:n]*pnds[:n]
    ptemp_n_inits = 1 - pinits[:n]

    for i in range(1,n+1):
        for j in range(i):
            temp_arr = np.array([1-_reinit_func(starts[k]-ends[j]) for k in range(j+1, i)])
            total_prenit = _reinit_func(starts[i]-ends[j])*np.product(temp_arr) # np.prod is 1.0 if temp_arr is empty
            pr[i] += pr[j]*ptemp_inits[j]*total_prenit
        pr[i] += np.product(ptemp_n_inits[:i])

    preach = pr[n]
    ptran = preach*pinits[n]*pelons[n]

    return ptran, preach


def wrapper1(orf_fasta, genome_fasta):
    """
    FUNCTION SIGNATURE
    ------------------
    param orf_fasta: str
        path to file containing the file containing all the orfs
    param genome_fasta: str
        path to file containing the file containing all the constructs
    return cds_tran: dict
        key -> CDS, value -> ptran
    return cds_reach: dict
        key -> CDS, value -> preach
    return maxOrfs: int
        maximum number of ORFs in a construct
    """
    (orf_pos, orf_dict_arranged, orf_fourth, pcomp_dict,
        genome_ids, input_cdss) = _generate_inputs(orf_fasta, genome_fasta)
    cds_tran = {}
    cds_reach = {}
    maxOrfs = 0
    for ids in genome_ids:
        orfs_list = np.array(orf_dict_arranged[ids])
        input_cdss = np.array(input_cdss)
        arr = [orf in input_cdss for orf in orfs_list]
        assert(sum(arr) == 1)
        n = arr.index(True)
        print(orfs_list[n], n)
        if n>maxOrfs:
            maxOrfs = n
        ptran_arr = np.zeros(n+1) # original, 1st deleted, 2nd deleted, ..., nth deleted
        preach_arr = np.zeros(n+1)
        starts = np.array([orf_pos[orf][0] for orf in orfs_list[:n+1]])
        ends = np.array([orf_pos[orf][1] for orf in orfs_list[:n+1]])
        pinits = np.array([_ppos_func(orf_fourth[orf]) for orf in orfs_list[:n+1]])
        pelons = np.array([pcomp_dict[orf] for orf in orfs_list[:n+1]])
        pnds = np.array([_nd_func(ends[i]-starts[i]) for i in range(n+1)])
        ptran, preach = find_ptran(n, pinits, pelons, pnds, starts, ends)
        ptran_arr[0] = ptran
        preach_arr[0] = preach

        # finding effect of ORFs on ptran
        for i in range(n):
            # will delete one orf and find ptran
            starts_del = np.append(starts[:i], starts[i+1:])
            ends_del = np.append(ends[:i], ends[i+1:])
            pinits_del = np.append(pinits[:i], pinits[i+1:])
            pelons_del = np.append(pelons[:i], pelons[i+1:])
            pnds_del = np.append(pnds[:i], pnds[i+1:])
            ptran_del, preach_del = find_ptran(n-1, pinits_del, pelons_del, pnds_del, starts_del, ends_del)
            ptran_arr[i+1] = ptran_del
            preach_arr[i+1] = preach_del

        cds_tran[orfs_list[n]] = ptran_arr
        cds_reach[orfs_list[n]] = preach_arr

    return cds_tran, cds_reach, maxOrfs

import sys
assert(len(sys.argv) == 4)
orf_fasta = sys.argv[1]
genome_fasta = sys.argv[2]
output_file = sys.argv[3]
assert(output_file[-3:]=="csv")

cds_tran, cds_reach, maxOrfs = wrapper1(orf_fasta, genome_fasta)

header = np.array(['CDS', 'Orig Ptran', 'Orig Preach'])
header = np.append(header, np.array([['Ptran - '+str(i+1)+' ORF', 'Preach - '+str(i+1)+' ORF'] for i in range(maxOrfs)]))

with open(output_file, 'w') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerow(header)
    for key, vals in cds_tran.iteritems():
        reach_vals = cds_reach[key]
        row = np.array([key])
        row = np.append(row, np.array([[vals[i], reach_vals[i]] for i in range(len(vals))]))
        csvwriter.writerow(row)

print("Completed!")
