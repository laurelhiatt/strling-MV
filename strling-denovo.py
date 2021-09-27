import argparse
import denovo

def get_args():
    """Incorporating argparse into the code for interchangeable arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--outliers", required=True,
    ### strling output goes here as input
        help="input outlier name")
    parser.add_argument("--ped", required=True,
        help="input ped file")
    ### ped file to sort trios
    parser.add_argument("--out",
        help="outputfile")
    ### out will just be the name of my output file... turn up
    parser.add_argument("--wiggle", default=0.1,
        help="establishes the ranges we are setting for the alleles (default: %(default)s)")
    ### cannot be greater than 1
    parser.add_argument("--minwig", default=10.0,
        help="minimum wiggle for small alleles")
    parser.add_argument("--depth", type=int, default=15,
        help="depth filter")
    parser.add_argument("--ampsize", type=int, default=150,
        help="amplification size filter")
        ### size of de novo expansion, or difference from kid to mom and dad allele sizes, is defaulted to 150bp
    return parser.parse_args()

def main():
    args = get_args()
    denovo.get_denovos(args)

if __name__ == "__main__":
	main()
