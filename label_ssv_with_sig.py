import argparse
import sys
from collections import defaultdict
from presig import presig as ps


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
            help="A MAF file (or MAF-like file)", required=True, dest="maf")
    parser.add_argument("-s", "--signatures",
            help="A Tidy Signature file of SBS signatures.", required=True, dest="tidysigs")
    parser.add_argument("-r", "--ref",
            help="A FASTA reference to retrieve flanking context from.", required=True, dest="ref")
    return parser.parse_args()


"""
Takes a tidy signature dictionary
and returns a new dictionary mapping from
a string (Context + "_" Change) -> most probable signature.

restrictToSigs is a list that can restrict to specific signatures (e.g., SBS1 and SBS5 ["SBS1", "SBS5"]
"""
def precompute_likely_sigs(tidy_sig_d, restrictToSigs=[]):
    feature_to_sig = defaultdict(str)

    inter_sig_d = defaultdict(str)
    if len(restrictToSigs) > 0:
        inter_sig_d = {k:tidy_sig_d[k] for k in tidy_sig_d if k in restrictToSigs}
    else:
        inter_sig_d = tidy_sig_d

    
    feat_sig_d = defaultdict(lambda : defaultdict(float))
    for sig in inter_sig_d:
        for change in inter_sig_d[sig]:
            for context in inter_sig_d[sig][change]:
                feat = ps.context_and_change_to_sbs_feature(context, change)
                feat_sig_d[feat][sig] = inter_sig_d[sig][change][context]
    
    for i in feat_sig_d:
        max_sig = max(feat_sig_d[i].keys(), key=(lambda k:feat_sig_d[i][k]))
        feature_to_sig[i] = max_sig
    return feature_to_sig

def calculate_probable_sig(sbs_feature, feature_to_sig):
    
    if sbs_feature in feature_to_sig:
        return feature_to_sig[sbs_feature]
    return None



if __name__ == "__main__":
    
    args = parse_args()

    ## Precompute the most likely signature for each feature
    tidy_sig_d = ps.parse_tidy_sig_file(args.tidysigs)[0]
    feature_to_sig = precompute_likely_sigs(tidy_sig_d, ["SBS1", "SBS5", "SBS13", "SBS2"])

    bc = ps.BasicContextualizer()
    bc.initialize_fasta(args.ref)


    header_d = {}
    sample_to_sig_to_count = defaultdict(lambda : defaultdict(int))
    ## Open MAF
    with open(args.maf, "r") as ifi:
        for line in ifi:
            ## For each line in the MAF file:
            tokens = line.strip().split("\t")
            if not "Chr" in line:
                ## Process only SNVs
                if tokens[header_d["Variant_Type"]] != "SNP":
                    continue
                chrom = tokens[header_d["Chromosome"]]
                start_pos = int(tokens[header_d["Start_position"]])
                ref_allele = tokens[header_d["Reference_Allele"]]
                alt_allele = tokens[header_d["Tumor_Seq_Allele"]]

                ## calculate an end_position if one is not provided
                endp = None
                if not "End_position" in header_d:
                    endp = ps.calculate_end_position(start_pos, ref_allele, alt_allele)
                else:
                    endp = tokens[header_d["End_position"]]

                ref_context_fiveprime, ref_context_threeprime = bc.get_reference_contexts(chrom, start_pos - 1, endp, 4)
                ## Produce an SBS feature for the SSV
                feat = ps.classify_SBS_feature(ref_allele, alt_allele, ref_context_fiveprime, ref_context_threeprime)
                ## Determine most likely signature of origin
                sample_to_sig_to_count[tokens[header_d["Patient_ID"]]][feature_to_sig[feat]] += 1

            else:
                ## Handle the header
                for i in range(0, len(tokens)):
                    header_d[tokens[i]] = i

    for i in sample_to_sig_to_count:
        for j in sample_to_sig_to_count[i]:
            print("\t".join([str(i), str(j), str(sample_to_sig_to_count[i][j])]))
