import sys
sys.path.append("..")
from denovo import *
import pytest
import argparse
import numpy as np

args = get_args(['--outliers', 'test.tsv', '--ped', 'test.ped'])
# additional args set with default values

@pytest.mark.parametrize("allele1, allele2, expected", [
    (100, 100, (100, 100)),  # check to see if allele baseline functions
    (0, 600, (0, 350)),      # check to see if large allele 2 is replaced
    (np.nan, 100, (100, 100)), # check to see if nan is replaced
    (100, np.nan, (100, 100)), # check that nan is replaced on the other side
    (np.nan, 400, (350, 350)), # check to see if nan is replaced at cutoff
    (400, np.nan, (350, 350)), # check to see if other nan is replaced at cutoff
    (400, 600, (350, 350))   # check to see if both large alleles are replaced
    ])

    # can't figure out how to get NaN to pytest work... code definitely works

def test_allele_check(allele1, allele2, expected):
    assert allele_check(allele1, allele2, args) == expected

@pytest.mark.parametrize("allele, expected", [
    (100, (90, 110.00000000000001)), # check to 0.1 easily, some rounding
    (0, (-10, 10)),                  # check to see how negative #s work
    (10, (0, 20)),                  # this and above check minwig, default 10
    (300, (270, 330))               # check to make sure 0.1 works with larger #
    ])

def test_wiggle(allele, expected):
    assert wiggle(allele, args) == expected

#def test_wiggle_error(allele):
#    with pytest.raises(ValueError):
#        wiggle(allele, args)
# No longer relevant here with args change of input

@pytest.mark.parametrize("allele1, allele2, expected", [
    (100, 100, ((90, 110.00000000000001), (90, 110.00000000000001))),
    # same alleles get same ranges
    (0, 100, ((-10, 10), (90, 110.00000000000001))),
    # different alelles get different ranges, with minwig included
    ])

def test_allele_range(allele1, allele2, expected):
    assert allele_range(allele1, allele2, args) == expected

@pytest.mark.parametrize("allele1, allele2, kidallele, expected", [
    (100, 100, 100, True), # kidallele matching parent allele returns True
    (np.nan, 100, 100, True), # nan allele doesn't change match
    (np.nan, 50, 45, True), # allele still matches if within range
    (500, 100, 1000, True), # large alleles match even if outside range
    (500, 100, 30, False), # if kid alelle is smaller, it is a deletion
    (30, 25, 50, False), # checking amplification
    (200, 300, 500, False), # check large allele amplification ***
    (0, np.nan, 3000, False), # large in kid vs 0 in parent
    ])

def test_check_range(allele1, allele2, kidallele, expected):
    assert check_range(allele1, allele2, kidallele, args) == expected


@pytest.mark.parametrize("mom_dict, dad_dict, kid_dict, expected", [
    ({'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150},
    {'allele1': 150, 'allele2': 150}, ('Full match', False)),
    # exact allele match across the board, should work

    ({'allele1': 3000, 'allele2': 150}, {'allele1': 150, 'allele2': 150},
    {'allele1': 150, 'allele2': 450}, ('Full match', False)),
    # match with large alleles

    ({'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150},
    {'allele1': 160, 'allele2': 140}, ('Full match', False)),
    # match within range, non-exact alleles

    ({'allele1': np.nan, 'allele2': np.nan}, {'allele1': 150, 'allele2': 150},
    {'allele1': 150, 'allele2': 150}, ('Missing alleles, ignore', False)),
    #ensure missing alleles are ignored

    ({'allele1': 0, 'allele2': np.nan}, {'allele1': np.nan, 'allele2': 0},
    {'allele1': 3000, 'allele2': 150}, ('Double MV, likely error', False)),
    # no kid alleles match a parent alllele, not an amp from DMV default

    ({'allele1': 20, 'allele2': 31.27}, {'allele1': 20.0, 'allele2': 88},
    {'allele1': 20, 'allele2': 36.26}, ('Full match', False)),
    # checking minwiggle and different allele order, kid2 to mom2

    ({'allele1': 0.0, 'allele2': 89.0}, {'allele1': 0.0, 'allele2': np.nan},
    {'allele1': 0.0, 'allele2': 80}, ('Full match', False)),
    # nan replaced  in dad and second allele matching with mom

    ({'allele1': 4000, 'allele2': 600}, {'allele1': 350, 'allele2': 555},
    {'allele1': 400, 'allele2': 1000}, ('Full match', False)),
    # large alleles over threshold should match

    ({'allele1': -1.0, 'allele2': np.nan}, {'allele1': -1.0, 'allele2': 90.92},
    {'allele1': -1.0, 'allele2': 117.48}, ('MV', False)),
    # close but no cigar, second allele outside of standard wiggle range

    ({'allele1': -1.0, 'allele2': np.nan}, {'allele1': -1.0, 'allele2': 90.92},
    {'allele1': -1.0, 'allele2': 300}, ('MV', True)),
    # Amplificatioin of this MV is true

    ({'allele1': 0, 'allele2': np.nan}, {'allele1': 0, 'allele2': 82.74},
    {'allele1': 0, 'allele2': np.nan}, ('Full match', False)),
    # ensuring allele check happens on kidalleledict
    ])

def test_full_allele_check(mom_dict, dad_dict, kid_dict, expected):
    assert full_allele_check(mom_dict, dad_dict, kid_dict, args) == expected
