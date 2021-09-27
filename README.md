# strling-MV

Necessary python modules:
argparse
pandas
numpy
peddy (must be in python version 3.7 for peddy)
math


Install dependencies using conda:
`conda env create environment.yml`  
`conda activate strling-denovo`

Utilizing STRling output as input, code now counts amplifications over bp threshold and determines Mendelian status (matched to parents, violation, etc.) at a certain depth threshold.

Command line might look something like:

python strling-MV.py --outliers STRs.tsv --ped file.ped --out output.tsv --wiggle 0.3 --minwig 10 --ampsize 100 --depth 12
