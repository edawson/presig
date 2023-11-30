import argparse
import sys
from collections import defaultdict
from pyfaidx import Fasta


class FastaCache:
    fasta = None
    seqcache = None
    max_cache_size = None

    def __init__(self):
        self.fasta = None
        self.seqcache = defaultdict(str)
        self.max_cache_size = 5
    
    
    def initialize_fasta(self, fa):
        self.fasta = Fasta(fa)
        self.seqcache = defaultdict(str)

    def get_seq(self, seqname):
        seqname = str(seqname)  
        if not seqname in self.seqcache:
            self.seqcache.clear()
            self.seqcache[seqname] = str(self.fasta[seqname])
        return self.seqcache[seqname]
    
    def get_seq_substr(self, seqname, start, end):
        if not seqname in self.seqcache:
            self.seqcache.clear()
            self.seqcache[seqname] = str(self.fasta[seqname])
        return self.seqcache[seqname][start:end]
    
    def get_reference_contexts(self, seqname, zero_based_start, zero_based_end, context_len = 25):
        ref_context_fiveprime =  self.get_seq_substr(seqname, zero_based_start - context_len, zero_based_start)
        ref_context_threeprime = self.get_seq_substr(seqname, zero_based_end + 1, zero_based_end + 1 + context_len)
        return ref_context_fiveprime, ref_context_threeprime

class BasicContextualizer:

    def __init__(self):
        self.fasta = None
    
    def initialize_fasta(self, fa):
        self.fasta = Fasta(fa)

    def get_seq(self, seqname):
        return str(self.fasta[seqname])
    
    def get_subseq(self, seqname, start, end):
        return str(self.fasta[seqname][start:end])
    
    def get_reference_contexts(self, seqname, zero_based_start, zero_based_end, context_len = 25):
        ref_context_fiveprime = str(self.fasta[seqname][zero_based_start - context_len:zero_based_start])
        ref_context_threeprime = str(self.fasta[seqname][zero_based_end+1:zero_based_end+1 + context_len])
        return ref_context_fiveprime, ref_context_threeprime

class CachedContextualizer:

    def __init__(self):
        self.fc = None
    
    def initialize_fasta(self, fa):
        self.fc = FastaCache()
        self.fc.initialize_fasta(fa)
    
    def get_seq(self, seqname):
        return self.fc.get_seq(seqname)
    
    def get_subseq(self, seqname, start, end):
        return self.fc.get_seq_substr(seqname, start, end)
    
    def get_reference_contexts(self, seqname, zero_based_start, zero_based_end, context_len = 25):
        return self.fc.get_reference_contexts(seqname, zero_based_start, zero_based_end, context_len)

def chop(s, ext):
    if s.endswith(ext):
        return s[:-len(ext)]
    return s
 
def basename(fi, ext = ".maf"):
    return chop(fi.split("/")[-1], ext)

def dirname(fi):
    return "/".join(fi.split("/")[:-1])

GLOBAL_MAX_MH_LEN = 5
GLOBAL_MAX_RPT_LEN = 5
GLOBAL_MAX_INDEL_LEN = 5

GLOBAL_LT_INDEL_WARNING = False
GLOBAL_UNKOWN_TYPE_WARNING = False

GLOBAL_PYRIMIDINE_LIST = ["C", "T"]

GLOBAL_SBS96_CONTEXTS = ["ACA","ACC","ACG","ACT","ACA","ACC",
                     "ACG","ACT", "ACA","ACC","ACG","ACT",
                     "ATA","ATC","ATG","ATT","ATA","ATC",
                     "ATG","ATT","ATA","ATC","ATG","ATT",
                     "CCA","CCC","CCG","CCT","CCA","CCC",
                     "CCG","CCT","CCA","CCC","CCG","CCT",
                     "CTA","CTC","CTG","CTT","CTA","CTC",
                     "CTG","CTT","CTA","CTC","CTG","CTT",
                     "GCA","GCC","GCG","GCT","GCA","GCC",
                     "GCG","GCT","GCA","GCC","GCG","GCT",
                     "GTA","GTC","GTG","GTT","GTA","GTC",
                     "GTG","GTT","GTA","GTC","GTG","GTT",
                     "TCA","TCC","TCG","TCT","TCA","TCC",
                     "TCG","TCT","TCA","TCC","TCG","TCT",
                     "TTA","TTC","TTG","TTT","TTA","TTC",
                     "TTG","TTT","TTA","TTC","TTG","TTT"]

GLOBAL_SBS96_CHANGES = ["C>A","C>G","C>T",
                    "T>A","T>C","T>G"]
        

GLOBAL_SBS96_FEATURES = [
"A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G",
"A[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","A[T>A]A","A[T>A]C",
"A[T>A]G","A[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","A[T>G]A",
"A[T>G]C","A[T>G]G","A[T>G]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T",
"C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","C[C>T]A","C[C>T]C","C[C>T]G",
"C[C>T]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","C[T>C]A","C[T>C]C",
"C[T>C]G","C[T>C]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[C>A]A",
"G[C>A]C","G[C>A]G","G[C>A]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T",
"G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","G[T>A]A","G[T>A]C","G[T>A]G",
"G[T>A]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","G[T>G]A","G[T>G]C",
"G[T>G]G","G[T>G]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","T[C>G]A",
"T[C>G]C","T[C>G]G","T[C>G]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
"T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","T[T>C]A","T[T>C]C","T[T>C]G",
"T[T>C]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"
]


GLOBAL_ID83_FEATURES = [
"1:Del:C:0","1:Del:C:1","1:Del:C:2","1:Del:C:3","1:Del:C:4","1:Del:C:5",
"1:Del:T:0","1:Del:T:1","1:Del:T:2","1:Del:T:3","1:Del:T:4","1:Del:T:5",

"1:Ins:C:0","1:Ins:C:1","1:Ins:C:2","1:Ins:C:3","1:Ins:C:4","1:Ins:C:5",
"1:Ins:T:0","1:Ins:T:1","1:Ins:T:2","1:Ins:T:3","1:Ins:T:4","1:Ins:T:5",

"2:Del:R:0","2:Del:R:1","2:Del:R:2","2:Del:R:3","2:Del:R:4","2:Del:R:5",
"3:Del:R:0","3:Del:R:1","3:Del:R:2","3:Del:R:3","3:Del:R:4","3:Del:R:5",
"4:Del:R:0","4:Del:R:1","4:Del:R:2","4:Del:R:3","4:Del:R:4","4:Del:R:5",
"5:Del:R:0","5:Del:R:1","5:Del:R:2","5:Del:R:3","5:Del:R:4","5:Del:R:5",

"2:Ins:R:0","2:Ins:R:1","2:Ins:R:2","2:Ins:R:3","2:Ins:R:4","2:Ins:R:5",
"3:Ins:R:0","3:Ins:R:1","3:Ins:R:2","3:Ins:R:3","3:Ins:R:4","3:Ins:R:5",
"4:Ins:R:0","4:Ins:R:1","4:Ins:R:2","4:Ins:R:3","4:Ins:R:4","4:Ins:R:5",
"5:Ins:R:0","5:Ins:R:1","5:Ins:R:2","5:Ins:R:3","5:Ins:R:4","5:Ins:R:5",

"2:Del:M:1","3:Del:M:1","3:Del:M:2","4:Del:M:1","4:Del:M:2","4:Del:M:3",
"5:Del:M:1","5:Del:M:2","5:Del:M:3","5:Del:M:4","5:Del:M:5"
]

GLOBAL_ID83_HASHSET = set(GLOBAL_ID83_FEATURES)
GLOBAL_SBS96_HASHSET = set(GLOBAL_SBS96_FEATURES)

"""
A debugging function.
Given a list of values with any type,
write those values to stderr separated by spaces.
"""
def write_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

"""
Creates a new header dictionary from a tab-separated header.
This can then be used to get the index of a given column name.
"""
def reheader(headerline):
    header = defaultdict(int)
    headertokens = headerline.strip("#").strip().split("\t")

    for i in range(0, len(headertokens)):
        header[headertokens[i]] = i

    return header 

"""
Reads a file with the following columns:
Signature Change  Context Amount

and generates a dictionary of signature:change:context:amount
(i.e., tidy_sig_d["SBS1]["C>T"]["ACG"]=0.3)
as well as
a dictionary mapping signature names to a vector of probabilities for each change/context

If isSampleData is passed as true, no validation is performed.
Validates that the number of dictionary entries is equal to nfeatures.
"""
def parse_tidy_sig_file(sig_file, isSampleData = False, nFeatures = 96):
    nested_d = lambda : defaultdict(nested_d)
    tidy_sig_d = nested_d()
    vec_d = defaultdict(list)

    header_d = None
    r_header_d = None

    with open(sig_file, "r") as ifi:
        for line in ifi:
            line = line.strip()
            tokens = line.split("\t")
            if not "Context" in line:
                tidy_sig_d[ tokens[header_d["Signature"]] ][ tokens[header_d["Change"]] ][ tokens[header_d["Context"]] ] = float(tokens[header_d["Amount"]])
                #tidy_sig_d[ tokens[header_d["Signature"]] ] = float(tokens[header_d["Amount"]])

            else:
                if header_d is not None:
                    write_err("ERROR: Two headers present. Exiting.")
                    exit(9)
                header_d = defaultdict(str)
                r_header_d = defaultdict(int)
                ind = 0
                for i in tokens:
                    header_d[i] = ind
                    r_header_d[ind] = i
                    ind = ind + 1

        return tidy_sig_d, vec_d
 

"""
Undoes alignment trimming of indels,
replacing bases represented with a "-" with
the correct base in the reference
and prepending the reference base to the non-dash allele.
"""
def untrim_indel(chrom, start, zb_start, end, zb_end, ref, alt, fa_ref):
    ## Get the REF allele
    real_ref = str(fa_ref[chrom][zb_start:zb_end])
    if alt == "-":
        ref = real_ref + ref
        alt = real_ref
        #start = int(start) - 1
    elif ref == "-":
        ref = real_ref
        alt = real_ref + alt
        #start = int(start) - 1
    else:
        write_err(["Invalid allele: ", ref, "->", alt, "."])
        raise Exception("Exiting.")
    return chrom, start, end, ref, alt

def context_and_change_to_sbs_feature(context, change):
    return join_sbs(context[0], change[0], change[-1], context[-1])

"""
Produces a trinucleotide context feature string
from a ref/alt allele and the 5' and 3' bases.
"""
def join_sbs(pre, ref, alt, post, quote = False):
    if quote:
        return "".join([ "\"", pre, "[", ref, ">", alt, "]", post, "\""])
    else:
        return "".join([pre, "[", ref, ">", alt, "]", post])

"""
Produces an indel feature string
given the set of indel feature attributes.
"""
def join_id(length, vtype, context, context_length, quote = False):
    if quote:
        return "".join(["\"", str(length), ":", str(vtype), ":", context, ":", str(context_length), "\""])
    return "".join([str(length), ":", str(vtype), ":", str(context), ":", str(context_length)])

def write_matlab_csv(d, features, filename):
    return
    

"""
Given a dictionary of features, a list of features
and an output filename
write a SigProfilerExtractor-compatible TSV file.
"""
def write_python_tsv(d, features, header, filename):
    with open(filename, "w") as ofi:
        ofi.write(header + "\n")
        for x in features:
            lab = x
            vals = [d[i][x] for i in sorted(d.keys())]
            line = "\t".join([lab, "\t".join(str(v) for v in vals)+ "\n"])
            ofi.write(line)

"""
Given a label for the feature column and a list of samples
creates the correct header
"""
def make_sample_header(feature_label, samples, sep = "\t", quote = False):
    return sep.join([feature_label, sep.join(sorted(samples))])


"""
Returns the length of an indel,
which is the difference in size between
its REF and ALT alleles. We strip off any dashes,
as they are used to indicate a single ref base usually
and unnecessarily make counting harder. I have never observed a MAF
with multipel dashes per line.
"""
def calculate_indel_length(ref, alt):
    ref = ref.strip("-")
    alt = alt.strip("-")
    return abs(len(ref) - len(alt)), len(ref), len(alt)

"""
Given a reference start position and a reference + alternate allele,
calculate the genomic end position of the variant.
"""
def calculate_end_position(start_position, ref, alt):
    var_len, ref_len, alt_len = calculate_indel_length(ref, alt)
    ## Set the var_len to 1 if we have an insertion OR
    ## a SNP
    if alt_len > ref_len:
        var_len = 1
    elif var_len == 0:
        var_len = alt_len
    return start_position + var_len - 1

"""
Creates a minimal representation of a variant that is compatible with the SigProfiler text format.
"""
def make_minimal_record(cancer_type, sample, assay, genome, variant_type, chrom, pos, end, ref, alt, mutation_type):
    vals = [cancer_type, sample, assay, genome, variant_type, chrom, pos, end, ref, alt, mutation_type]
    return "\t".join(vals)

"""
Reverse complements a sequence
"""
GLOBAL_REVCOMP_D = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def reverse_complement(seq):
    return "".join([GLOBAL_REVCOMP_D[i] for i in seq[::-1]])


"""
Strand complements a sequence (i.e., converts
it to the canonical T/C strand).
"""
GLOBAL_STRANDCOMP_D = {'A':'T','C':'C','G':'C','T':'T','N':'N'}
def strand_complement(seq, reverse= False):
    return "".join([GLOBAL_STRANDCOMP_D[i] for i in seq]) if not reverse else "".join([GLOBAL_STRANDCOMP_D[i] for i in seq[::-1]])


"""
Handles command line argurment parsing with argparse
"""
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--maf", required=True, dest="maf", help="MAF file from which to extract matrix.")
    parser.add_argument("-i", "--indels", action="store_true", dest="indels", help="Extract indel matrix, as well as SBS", default=False)
    parser.add_argument("-s", "--sigprofiler", help = "Output in sigprofiler format", action = "store_true")
    parser.add_argument("-u", "--untrim", help = "Untrim indels before counting features or producing output.", action = "store_true", dest="untrim")
    parser.add_argument("-R", "--addref", help = "Attach a reference base to the indel call, but don't untrim.", action = "store_true", dest="addref")
    parser.add_argument("-e", "--exome", dest="exome", help="Exome data - restrict genome to exome regions", default=False)
    parser.add_argument("-p", "--project", dest="project", default="PROJECT", help="Project name for output.")
    parser.add_argument("-d", "--directory", dest="directory", default="input", help="Input/Output directory")
    parser.add_argument("-f", "--fasta", required=True, dest="ref", help="A fasta file reference.")
    parser.add_argument("-C", "--custom-id", dest="custom_id", help="Column name of custom, preferred ID field", default=None, type = str)
    parser.add_argument("-L", "--log", dest="log", action="store_true", help="Log each mutation, its feature and its context to a logfile.")

    return parser.parse_args()

"""
Detects microhomology on flanking portions of variants
"""
def detect_microhomology(ref_allele, alt_allele,
                            reflen, altlen,
                            ref_context_fiveprime, ref_context_threeprime):
    mh = False
    mh_len = 0
    head_mh_seq = ""
    tail_mh_seq = ""
    head_mh_valid = True
    tail_mh_valid = True
    for i in range(1, 1 + GLOBAL_MAX_MH_LEN):
        if reflen > altlen:
            head_mh_seq = ref_allele[0:i]
            tail_mh_seq = ref_allele[len(ref_allele) - i : len(ref_allele)]
        else:
            head_mh_seq = alt_allele[0:i]
            tail_mh_seq = alt_allele[len(alt_allele) - i : len(alt_allele)]
        fiveprime_mh_seq = ref_context_fiveprime[len(ref_context_fiveprime) - i : len(ref_context_fiveprime)]
        threeprime_mh_seq = ref_context_threeprime[0:i]
        if tail_mh_seq == fiveprime_mh_seq:
            # write_err("5\' symm:", tail_mh_seq, ";", "context:",
            #     ref_context_fiveprime,
            #     "ref:", ref_allele,
            #     "alt:", alt_allele)
            mlen = len(tail_mh_seq)
            mh_len = max(mh_len, mlen)
        if head_mh_seq == threeprime_mh_seq:
            mlen = len(head_mh_seq)
            mh_len = max(mh_len, mlen)
    if mh_len > 0:
        mh = True

    return mh, mh_len

"""
Handles repeat counts for any number of base pairs in both directions
using a single loop.
Takes the ref and alt alleles, their lengths, and the fiveprime/threeprime ref contexts.
Returns a boolean indicating whether the sequence is an INDEL of a repeat unit
and the "len", which is really the number of repeat units surrounding the variant.
"""
def detect_repeat(ref_allele, alt_allele,
                    ref_len, alt_len,
                    ref_context_fiveprime, ref_context_threeprime):
    
    global GLOBAL_LT_INDEL_WARNING
    rpt = False
    rpt_len = 0
    ## Handle repeat counts in either direction,
    ## though most VCF / MAF files will be left-aligned and trimmed.
    ## First check repeats 3' of our SSV
    seq = ""

    ## This line handles both trimmed and untrimmed alleles
    ## by removing dashes and any remaining bases.
    seq = ref_allele[len(alt_allele.strip("-")):] if ref_len > alt_len else alt_allele[len(ref_allele.strip("-")):]

    seq_start = 0
    seq_end = len(seq)
    fiveprime_symmetry_valid = True
    for i in range(0, GLOBAL_MAX_RPT_LEN):
        ref_seq = ref_context_threeprime[seq_start:seq_end]
        fiveprime_ref_seq = ref_context_fiveprime[len(ref_context_fiveprime) - seq_end: len(ref_context_fiveprime) - seq_start]
        #write_err(fiveprime_ref_seq, seq, ref_seq)                     
        if (fiveprime_symmetry_valid and fiveprime_ref_seq == seq):
            rpt = True
            if not GLOBAL_LT_INDEL_WARNING:
                write_err("Warning: INDELS don't appear to be left-aligned.", ref_context_fiveprime, ref_allele, alt_allele, ref_context_threeprime)
                GLOBAL_LT_INDEL_WARNING = True
            rpt_len = rpt_len + 1
            if rpt_len == GLOBAL_MAX_RPT_LEN:
                break
        else:
            fiveprime_symmetry_valid = False

        #write_err(start_pos, seq, ref_seq)
        if ref_seq == seq:
            rpt = True
            rpt_len = rpt_len + 1
            if rpt_len == GLOBAL_MAX_RPT_LEN:
                break
        else:
            break
        seq_start += len(seq)
        seq_end += len(seq)                           

    return rpt, rpt_len

def classify_SBS_feature(ref_allele, alt_allele, ref_context_fiveprime, ref_context_threeprime):
    sbs_feature = None

    fiveprime_base = ref_context_fiveprime[-1] if ref_allele in GLOBAL_PYRIMIDINE_LIST else reverse_complement(ref_context_threeprime[0])
    threeprime_base = ref_context_threeprime[0] if ref_allele in GLOBAL_PYRIMIDINE_LIST else reverse_complement(ref_context_fiveprime[-1])

    if ref_allele not in GLOBAL_PYRIMIDINE_LIST:
        alt_allele = reverse_complement(alt_allele)
        ref_allele = reverse_complement(ref_allele)
    sbs_feature = join_sbs(fiveprime_base, strand_complement(ref_allele), alt_allele, threeprime_base)
    assert sbs_feature in GLOBAL_SBS96_HASHSET
    
    return sbs_feature


"""
Determines what feature a given variant represents (e.g. SBS A[C>T]G, or one
of the ID83 features, etc).
Takes 
- a variant (as a line in a MAF file)
- a contextualizer, which is just a wrapper around how to get
the sequence context of a variant and
- a header for determining which fields to grab.
It then returns a feature and a feature type in ["ID", "SBS", "DBS"]
"""
def maf_line_to_feature(line,
                        contextualizer,
                        header_d,
                        sbs_d,
                        id_d,
                        logfi = None,
                        sigprofiler = False,
                        args = None):
    ## The reported feature, e.g. A[C->A]G or 2:Del:M1
    feature = None
    ## One of indel, SBS, DBS, or SV
    feature_type = None

    line = line.strip().split("\t")
    vtype = tokens[header_d["Variant_Type"]]
    chrom = tokens[header_d["Chromosome"]]
    strand = tokens[header_d["Strand"]]
    sample = tokens[header_d[id_field]]
    start_pos = tokens[header_d["Start_position"]]
    zero_based_start = int(start_pos) - 1;
    end_pos = tokens[header_d["End_position"]]
    zero_based_end = int(end_pos) - 1;

    ref_allele = str(tokens[header_d["Reference_Allele"]]).upper()
    ref_len = len(ref_allele)
    alt_allele = str(tokens[header_d["Tumor_Seq_Allele2"]]).upper()

    #fasta_allele = str(ref[chrom][zero_based_start:zero_based_start+ref_len])
    assert(zero_based_start < zero_based_end + ref_len)
    fasta_allele = contextualizer.get_subseq(chrom, zero_based_start, zero_based_end + ref_len)

    strand = 1

    ref_context_fiveprime, ref_context_threeprime = contextualizer.get_reference_contexts(chrom, zero_based_start, zero_based_end, 25)
    #write_err(ref_context_fiveprime, ref_allele, fasta_allele, alt_allele, ref_context_threeprime)

    if vtype == "SNP" or vtype == "SNV":
        assert fasta_allele == ref_allele
        feature = classify_SBS_feature(ref_allele, alt_allele, ref_context_fiveprime, ref_context_threeprime)
        feature_type = "SBS"
        sbs_d[sample][feature] += 1

    elif vtype == "DEL" or vtype == "INS":
        ## Get the length of the variant
        vlen, reflen, altlen = calculate_indel_length(ref_allele , alt_allele)

        feat_type = "Del" if  reflen > altlen else "Ins"
        feat_len = min(vlen, GLOBAL_MAX_INDEL_LEN)
        
        feat_context = "NA"
        feat_context_len = 0
        rpt_len = 0
        rpt = False
        mh_len = 0
        mh = False

        ## If reflen > altlen (i.e., the variant is a deletion),
        ## the feature is the single base from the REF allele which remains
        ## after stripping off len(ALT allele) bases. 
        ## This base must then be strand-complemented (i.e., must be T or C)
        ## The old way of doing this was like so. This method is slower per variant than the
        ## ternary operator.
        # if vlen == 1 and feat_type == "Del":
        #     feat_context = strandcomp(ref_allele[len(alt_allele.strip("-")):])
        # elif vlen == 1 and feat_type == "Ins":
        #     feat_context = strandcomp(alt_allele[len(ref_allele.strip("-")):])
        rpt, rpt_len = detect_repeat(ref_allele, alt_allele,
                                        reflen, altlen,
                                        ref_context_fiveprime, ref_context_threeprime)
        mh, mh_len = detect_microhomology(ref_allele, alt_allele,
                                        reflen, altlen,
                                        ref_context_fiveprime, ref_context_threeprime)

        if vlen == 1:
            feat_context = strand_complement(ref_allele[len(alt_allele.strip("-")):])  if \
                reflen > altlen else \
                strand_complement(alt_allele[len(ref_allele.strip("-")):])
            feat_context_len = rpt_len
        elif rpt:
            feat_context = "R"
            feat_context_len = rpt_len
        elif mh and vtype == "DEL":
            feat_context = "M"
            feat_context_len = mh_len
        else:
            feat_context = "R"
            feat_context_len = 0


        feature = join_id(feat_len, feat_type, feat_context, feat_context_len)
        feature_type = "ID"
        sample_ID83_d[sample][feature] += 1
    else:
        global GLOBAL_UNKOWN_TYPE_WARNING
        if not GLOBAL_UNKOWN_TYPE_WARNING:
            write_err("Invalid variant type:", vtype)
            GLOBAL_UNKOWN_TYPE_WARNING = True
        return "", "", vtype

    if feature not in GLOBAL_SBS96_HASHSET and feature not in GLOBAL_ID83_HASHSET:
        write_err("ERROR: invalid feature", feature,
        sample, chrom, start_pos, end_pos, vtype, ref_context_fiveprime,  ref_allele, alt_allele, ref_context_threeprime)

    if args.sigprofiler:
        print(make_minimal_record(args.project, sample, "WGS", "GRCh37", vtype, chrom, start_pos, end_pos, ref_allele, alt_allele, "SOMATIC"))

    if logfi is not None:
        logfi.write("\t".join([sample, chrom, feature_type, feature, start_pos, end_pos, vtype, ref_allele, alt_allele, ref_context_fiveprime, ref_context_threeprime]) + "\n")


    return feature, feature_type, vtype

def main():

    ## Holds the sample->mutational counts vectors
    ## Each Key is a sample name
    sample_SBS96_d = defaultdict(lambda : defaultdict(int))
    sample_ID83_d = defaultdict(lambda : defaultdict(int))
    ## Holds the count of total variants processed (i.e., SNP: N, INS: N, DNP: N)
    total_var_count_d = defaultdict(int)


    header_d = {
    "Chromosome" : 4,
    "Start_position" : 5,
    "End_position" : 6,
    "Strand" : 7,
    "Variant_Type" : 9,
    "Reference_Allele" : 10,
    "Tumor_Seq_Allele1" : 11,
    "Tumor_Seq_Allele2" : 12,
    "Tumor_Sample_Barcode": 15,
    "Matched_Norm_Sample_Barcode" : 16
    }

    args = parse_args()

    header_line = None
    ref = None
    id_field = None

    if args.custom_id is not None:
        id_field = args.custom_id
    else:
        id_field = "Tumor_Sample_Barcode"

    logname = "logfile.txt"
    logfi = None
    if args.log:
        logfi = open("logfile.txt", "w")
    
    bc = BasicContextualizer()
    bc.initialize_fasta(args.ref)

    with open(args.maf, "r") as mfi:
        for line in mfi:
            line = line.strip()
            tokens = line.split("\t")
            if "Hugo_Symbol" not in line:
                feature, feature_type, vtype = maf_line_to_feature(line, bc,
                header_d, sample_SBS96_d,
                sample_ID83_d, logfi,
                args.sigprofiler, args)  
                total_var_count_d[vtype] += 1       
                
                # if vtype != "SNP" and vtype != "INS" and vtype != "DEL":
                #     write_err(["Invalid mutation type", vtype])

            else:
                header_line = line
                if (args.custom_id):
                    header_d = reheader(header_line)
                    continue
    
    write_err("Processed the following number of variants:")
    for t in total_var_count_d:
        write_err(t, ":", total_var_count_d[t])


    ## Write our SBS96 counts
    write_python_tsv(sample_SBS96_d, GLOBAL_SBS96_FEATURES, make_sample_header("MutationType", sample_SBS96_d.keys()), chop(basename(args.maf), ".tsv") + ".SBS96.tsv")
    ## Write our ID83 counts
    write_python_tsv(sample_ID83_d, GLOBAL_ID83_FEATURES, make_sample_header("MutationType", sample_ID83_d.keys()), chop(basename(args.maf), ".tsv") + ".ID83.tsv")
    if logfi is not None:
        logfi.close()
    




    

if __name__ == "__main__":
    main()
