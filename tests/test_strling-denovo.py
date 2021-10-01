import sys
sys.path.append("..")
from denovo import *
import pytest
import argparse

@pytest.mark.parametrize("allele1, allele2, expected", [
    (100, 100, (100, 100)),
])
def test_allele_check(allele1, allele2, expected):
    assert allele_check(allele1, allele2) == expected

@pytest.mark.skip
@pytest.mark.parametrize("allele, expected", [
    (100, (90, 110.00000000000001)),
    (0, (-1, 1)),
    (10, (5, 15)),
])

@pytest.mark.skip
def test_wiggle(allele, expected):
    assert wiggle(allele) == expected

# wiggle should raise an error when wgl is not between 0 and 1
@pytest.mark.skip
@pytest.mark.parametrize("allele", [
    (100, (90, 110)),
])

@pytest.mark.skip
def test_wiggle_error(allele):
    with pytest.raises(ValueError):
        wiggle(allele, wgl, minwig)

@pytest.mark.skip
@pytest.mark.parametrize("allele1, allele2, expected", [
    (100, 100, (90, 110.00000000000001), (90, 110.00000000000001)),
])

@pytest.mark.skip
def test_allele_range(allele1, allele2, expected):
    assert allele_range(allele1, allele2) == expected

@pytest.mark.skip
@pytest.mark.parametrize("allele1, allele2, expected", [
    (100, 0.1, 1, (90, 110.00000000000001)),
    (0, 0.1, 1, (-1, 1)),
    (10, 0.1, 5, (5, 15)),
])

@pytest.mark.skip
def test_get_allele_ranges(allele1, allele2, expected):
    assert get_allele_ranges(allele1, allele2) == expected

@pytest.mark.skip
@pytest.mark.parametrize("allele1, allele2, kidallele, expected", [
    (100, 0.1, 1, (90, 110.00000000000001)),
    (0, 0.1, 1, (-1, 1)),
    (10, 0.1, 5, (5, 15)),
])

@pytest.mark.skip
def test_check_range(allele1, allele2, kidallele, expected):
    assert check_allele_range(allele1, allele2, kidallele) == expected

@pytest.mark.skip
@pytest.mark.parametrize("mom_dict, dad_dict, kid_dict, expected", [(10, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, 'Full match'),
(10, {'allele1': 3000, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, 'Full match'), (100, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, 'Full match'),
(5, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 160, 'allele2': 150}, 'MV'), (10, {'allele1': 'NaN', 'allele2': 'NaN'}, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, 'Missing'),
(10, {'allele1': 3000, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 3000, 'allele2': 150}, 'MV')])

@pytest.mark.skip
def full_allele_check(mom_dict, dad_dict, kid_dict, expected):
    assert full_allele_check(momalleledict, dadalleledict, kidalleledict) == expected
