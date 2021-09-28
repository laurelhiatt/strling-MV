import argparse
import denovo

def get_args():
    """Incorporating argparse into the code for interchangeable arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--outliers", required=True,
    # strling output goes here as input
        help="input outlier name")
    parser.add_argument("--ped", required=True,
        help="input ped file")
    # ped file to sort trios
    parser.add_argument("--out",
        help="outputfile")
    # out will just be the name of the output file... turn up
    parser.add_argument("--wiggle", default=0.1,
        help="establishes range for alleles (default:%(default)s)")
    # value between 0 and 1
    parser.add_argument("--minwig", default=10.0,
        help="minimum wiggle for small alleles (default: %(default)s)")
    parser.add_argument("--depth", type=int, default=15,
        help="depth filter (default: %(default)s)")
    parser.add_argument("--ampsize", type=int, default=150,
        help="amplification size filter (default: %(default)s)")
    parser.add_argument("--allelecutoff", type=float, default=350.0,
        help="cutoff for max allele size (default: %(default)s)")
        # size of de novo expansion, or difference from kid to mom/dad alleles
    return parser.parse_args()

def main():
    """Culminating function to establish amplifications and Mendelian status"""
    args = get_args()
    denovo.get_denovos(args)

if __name__ == "__main__":
	main()
