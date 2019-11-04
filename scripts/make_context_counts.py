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
        self.fasta = pyfaidx.FASTA(fa)
        self.seqcache.clear()

    def __init__(self, fa):
        self.fasta = None
        self.seqcache = defaultdict(str)
        self.max_cache_size = 5
        initialize_fasta(fa)

    def get_seq(self, seqname):
        if len(seqcache) >= max_cache_size:
            self.seqcache.clear()
        if not seqname in seqcache:
            seqcache[str(seqname)] = str(self.fasta[seqname])
    
    def get_seq_substr(self, seqname, start, end):
        return self.get_seq(seqname)[start:end]


def chop(s, ext):
    if s.endswith(ext):
        return s[:-len(ext)]
    return s
 
def basename(fi, ext = ".maf"):
    return chop(fi.split("/")[-1], ext)

def dirname(fi):
    return "/".join(fi.split("/")[:-1])

        

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

# "2:Del:R:0","2:Del:R:1","2:Del:R:2","2:Del:R:3","2:Del:R:4","2:Del:R:5", "2:Del:R:6",
# "3:Del:R:0","3:Del:R:1","3:Del:R:2","3:Del:R:3","3:Del:R:4","3:Del:R:5", "3:Del:R:6",
# "4:Del:R:0","4:Del:R:1","4:Del:R:2","4:Del:R:3","4:Del:R:4","4:Del:R:5", "4:Del:R:6",
# "5:Del:R:0","5:Del:R:1","5:Del:R:2","5:Del:R:3","5:Del:R:4","5:Del:R:5", "5:Del:R:6",


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
def strand_complement(seq):
    return "".join([GLOBAL_STRANDCOMP_D[i] for i in seq])



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

    return parser.parse_args()


if __name__ == "__main__":

    ## Cache oft-used functions locally for performance.
    jsbs = join_sbs
    jid = join_id
    c_indel_len = calculate_indel_length
    strandcomp = strand_complement
    revcomp = reverse_complement


    ## Holds the sample->mutational counts vectors
    ## Each Key is a sample name
    sample_SBS96_d = defaultdict(lambda : defaultdict(int))
    sample_ID83_d = defaultdict(lambda : defaultdict(int))


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
    if args.ref is not None:
        ref = Fasta(args.ref)
    if args.custom_id is not None:
        id_field = args.custom_id
    else:
        id_field = "Tumor_Sample_Barcode"

    logfi = open("logfile.txt", "w")
    

    with open(args.maf, "r") as mfi:
        for line in mfi:
            line = line.strip()
            tokens = line.split("\t")
            if "Hugo_Symbol" not in line:
                vtype = tokens[header_d["Variant_Type"]]
                chrom = tokens[header_d["Chromosome"]]
                sample = tokens[header_d[id_field]]
                start_pos = tokens[header_d["Start_position"]]
                zero_based_start = int(start_pos) - 1;
                end_pos = tokens[header_d["End_position"]]
                zero_based_end = int(end_pos) - 1;

                ref_allele = str(tokens[header_d["Reference_Allele"]]).upper()
                alt_allele = str(tokens[header_d["Tumor_Seq_Allele2"]]).upper()

                ref_len = len(ref_allele)
                fasta_allele = str(ref[chrom][zero_based_start:zero_based_start+ref_len])

                strand = 1

                sbs_feature = ""
                indel_feature = ""

                
                
                ref_context_fiveprime = str(ref[chrom][zero_based_start - 25:zero_based_start])
                ref_context_threeprime = str(ref[chrom][zero_based_end+1:zero_based_end+1 + 25])
                #write_err(ref_context_fiveprime, ref_allele, fasta_allele, alt_allele, ref_context_threeprime)

                if vtype == "SNP":
                    assert fasta_allele == ref_allele
                    if strandcomp(ref_allele) == alt_allele:
                        alt_allele = revcomp(alt_allele)
                    sbs_feature = jsbs(ref_context_fiveprime[-1], strandcomp(ref_allele), alt_allele, ref_context_threeprime[0])
                    #write_err(sbs_feature)
                    if args.sigprofiler:
                        ## make_minimal_record(cancer_type, sample, assay, genome, variant_type, chrom, pos, ref, alt, mutation_type = "Somatic"):
                        print(make_minimal_record(args.project, sample, "WGS", "GRCh37", vtype, chrom, start_pos, end_pos, ref_allele, alt_allele, "SOMATIC"))
                    else:
                        sample_SBS96_d[sample][sbs_feature] += 1
                    
                elif vtype == "DEL" or vtype == "INS":

                    ## Get the length of the variant
                    vlen, reflen, altlen = c_indel_len(ref_allele , alt_allele)
                    feat_type = "NA"
                    feat_len = min(vlen, 5)
                    feat_context = "NA"
                    feat_context_len = 0

                    max_rpt_len  = 5
                    rpt_len = 0
                    rpt = False

                    max_mh_len = 5
                    mh_len = 0
                    mh = False


                    ## Get the variant type
                    ## These checks are currently overly thorough and use unreachable conditions,
                    ## given that we require vtype == "DEL|INS" above. TODO
                    vlower = vtype.lower()
                    if reflen > altlen and vlower == "del" or vlower == "deletion":
                        feat_type = "Del"
                    elif reflen < altlen and vlower == "ins" or vlower == "insertion":
                        feat_type = "Ins"
                    else:
                        write_err("Undefined variant type", vtype, chrom, start_pos, sample)
                        exit(9)

                    ## TODO: handle strandedness



                    ## Handle repeat counts in either direction,
                    ## though most VCF / MAF files will be left-aligned and trimmed.
                    ## First check repeats 3' of our SSV
                    seq = ""
                    if vtype == "DEL":
                        ## This line handles both trimmed and untrimmed alleles
                        ## by removing dashes and any remaining bases.
                        seq = ref_allele[len(alt_allele.strip("-")):]
                    else:
                        seq = alt_allele[len(ref_allele.strip("-")):]
                    seq_start = 0
                    seq_end = len(seq)
                    fiveprime_symmetry_valid = True
                    for i in range(0, max_rpt_len):
                        ref_seq = ref_context_threeprime[seq_start:seq_end]
                        fiveprime_ref_seq = ref_context_fiveprime[len(ref_context_fiveprime) - seq_end: len(ref_context_fiveprime) - seq_start]
                        #write_err(fiveprime_ref_seq, seq, ref_seq)                     
                        if (fiveprime_symmetry_valid and fiveprime_ref_seq == seq):
                            rpt = True
                            write_err("Warning: INDELS don't appear to be trimmed.", vtype, ref_context_fiveprime, ref_allele, alt_allele, ref_context_threeprime)
                            rpt_len = rpt_len + 1
                            if rpt_len == max_rpt_len:
                                break
                        else:
                            fiveprime_symmetry_valid = False

                        #write_err(start_pos, seq, ref_seq)
                        if ref_seq == seq:
                            rpt = True
                            rpt_len = rpt_len + 1
                            if rpt_len == max_rpt_len:
                                break
                        else:
                            break
                        seq_start += len(seq)
                        seq_end += len(seq)
                    #write_err(start_pos, ref_seq, rpt, mh)
                    if rpt and vlen > 1:
                        feat_context = "R"
                        feat_context_len = rpt_len
                    if not mh:
                        feat_context_len = rpt_len
                        if vlen > 1:
                            feat_context = "R"
                    if (vtype == "INS" and feat_context_len >=5) or (vtype == "DEL" and vlen == 1 and feat_context_len >= 5):
                        feat_context_len = 5                            


                    
                    ## Handle the base feature of 1-base deletions / insertions
                    if vlen == 1 and feat_type == "Del":
                        feat_context = strandcomp(ref_allele[len(alt_allele.strip("-")):])
                        assert len(feat_context) == 1
                    elif vlen == 1 and feat_type == "Ins":
                        feat_context = strandcomp(alt_allele[len(ref_allele.strip("-")):])
                        #write_err(ref_allele, alt_allele, feat_context, len(feat_context))
                        assert len(feat_context) == 1

                    ## Handle microhomology
                    head_mh_seq = ""
                    tail_mh_seq = ""
                    head_mh_valid = True
                    tail_mh_valid = True
                    if vlen > 1 and vtype == "DEL" and not rpt:
                        for i in range(1, 1+max_mh_len):
                            if vtype == "DEL":
                                head_mh_seq = ref_allele[0:i]
                                tail_mh_seq = ref_allele[len(ref_allele) - i : len(ref_allele)]
                            else:
                                head_mh_seq = alt_allele[0:i]
                                tail_mh_seq = alt_allele[len(alt_allele) - i : len(alt_allele)]
                            fiveprime_mh_seq = ref_context_fiveprime[len(ref_context_fiveprime) - i : len(ref_context_fiveprime)]
                            threeprime_mh_seq = ref_context_threeprime[0:i]
                            if tail_mh_seq == fiveprime_mh_seq:
                                write_err("5\' symm:", tail_mh_seq, ";", "context:",
                                 ref_context_fiveprime,
                                  "ref:", ref_allele,
                                  "alt:", alt_allele)
                                mlen = len(tail_mh_seq)
                                mh_len = max(mh_len, mlen)
                            if head_mh_seq == threeprime_mh_seq:
                                mlen = len(head_mh_seq)
                                mh_len = max(mh_len, mlen)
                        if mh_len > 0:
                            feat_context = "M"
                            feat_context_len = mh_len
                            mh = True

                    
                    
                    indel_feature = jid(feat_len, feat_type, feat_context, feat_context_len)
                    #write_err(indel_feature)

                    if args.sigprofiler:
                        pass
                    sample_ID83_d[sample][indel_feature] += 1
                    

                else:
                    write_err(["Invalid mutation type", vtype])
                if indel_feature not in GLOBAL_ID83_HASHSET and sbs_feature not in GLOBAL_SBS96_HASHSET:
                    write_err("ERROR: invalid feature", sbs_feature, indel_feature,
                    sample, chrom, start_pos, end_pos, vtype, ref_allele, alt_allele)
                logfi.write("\t".join([sample, chrom, sbs_feature, indel_feature, start_pos, end_pos, vtype, ref_allele, alt_allele, ref_context_fiveprime, ref_context_threeprime]) + "\n")

            else:
                header_line = line
                if (args.custom_id):
                    header_d = reheader(header_line)
                    continue


    ## Write our SBS96 counts
    write_python_tsv(sample_SBS96_d, GLOBAL_SBS96_FEATURES, make_sample_header("MutationType", sample_SBS96_d.keys()), chop(basename(args.maf), ".tsv") + ".SBS96.tsv")
    ## Write our ID83 counts
    write_python_tsv(sample_ID83_d, GLOBAL_ID83_FEATURES, make_sample_header("MutationType", sample_ID83_d.keys()), chop(basename(args.maf), ".tsv") + ".ID83.tsv")
    logfi.close()
    




    
