import argparse
from os.path import dirname
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", help="MAF file from which to extract matrix.", required=True)
    parser.add_argument("-e", "--exome", dest="exome", help="Exome data - restrict genome to exome regions", default=False)
    parser.add_argument("-p", "--project", dest="project", default="PROJECT", help="Project name for output.")
    parser.add_argument("-d", "--directory", dest="directory", default="input", help="Input/Output directory")
    parser.add_argument("-P", "--plot", dest="plot", action="store_true", help="Output plots of input data.")

    return parser.parse_args()


if __name__ == "__main__":
    
    args = parse_args()

    if args.directory is None:
        args.directory = dirname(args.input)

    matrices = matGen.SigProfilerMatrixGeneratorFunc(args.project, "GRCh37", dirname(args.input), plot=args.plot, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=True, cushion=100)
