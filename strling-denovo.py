import sys
import argparse
import denovo

def main(commandlineargs):
    """Culminating function to establish amplifications and Mendelian status"""
    args = get_args(commandlineargs)
    denovo.get_denovos(args)

if __name__ == "__main__":
	main(sys.argv[1:])
