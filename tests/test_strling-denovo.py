import sys
sys.path.append("..")
from denovo import *
import pytest

@pytest.mark.parametrize("allele1, allele2, expected", [
    (100, 0.1, 1, (90, 110.00000000000001)),
    (0, 0.1, 1, (-1, 1)),
    (10, 0.1, 5, (5, 15)),
])
def test_allele_check(allele1, allele2, expected):
    assert allele_check(allele1, allele2) == expected

@pytest.mark.parametrize("allele, expected", [
    (100, 0.1, 1, (90, 110.00000000000001)),
    (0, 0.1, 1, (-1, 1)),
    (10, 0.1, 5, (5, 15)),
])
def test_wiggle(allele, expected):
    assert wiggle(allele) == expected

# wiggle should raise an error when wgl is not between 0 and 1
@pytest.mark.parametrize("allele", [
    (100, 10, 1),
    (100, -10, 1),
])
def test_wiggle_error(allele):
    with pytest.raises(ValueError):
        wiggle(allele, wgl, minwig)

@pytest.mark.parametrize("allele1, allele2, expected", [
    (100, 0.1, 1, (90, 110.00000000000001)),
    (0, 0.1, 1, (-1, 1)),
    (10, 0.1, 5, (5, 15)),
])

def test_allele_range(allele1, allele2, expected):
    assert allele_range(allele1, allele2) == expected

@pytest.mark.parametrize("allele1, allele2, expected", [
    (100, 0.1, 1, (90, 110.00000000000001)),
    (0, 0.1, 1, (-1, 1)),
    (10, 0.1, 5, (5, 15)),
])

def test_get_allele_ranges(allele1, allele2, expected):
    assert get_allele_ranges(allele1, allele2) == expected

@pytest.mark.parametrize("allele1, allele2, kidallele, expected", [
    (100, 0.1, 1, (90, 110.00000000000001)),
    (0, 0.1, 1, (-1, 1)),
    (10, 0.1, 5, (5, 15)),
])

def test_check_range(allele1, allele2, kidallele, expected):
    assert check_allele_range(allele1, allele2, kidallele) == expected


@pytest.mark.parametrize("mom_dict, dad_dict, kid_dict, expected", [(10, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, 'Full match'),
(10, {'allele1': 3000, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, 'Full match'), (100, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, 'Full match'),
(5, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 160, 'allele2': 150}, 'MV'), (10, {'allele1': 'NaN', 'allele2': 'NaN'}, {'allele1': 150, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, 'Missing'),
(10, {'allele1': 3000, 'allele2': 150}, {'allele1': 150, 'allele2': 150}, {'allele1': 3000, 'allele2': 150}, 'MV')])

def full_allele_check(mom_dict, dad_dict, kid_dict, expected):
    assert full_allele_check(momalleledict, dadalleledict, kidalleledict) == expected
# to do
#@pytest.mark.parametrize("allele, wgl, minwig, expected", [
#    (10, 0.1, (1,500), (1,500) ),
#])
#def test_check_range(minwig, proportion, allele1, allele2, kidallele, expected):
#    assert check_range(minwig, proportion, allele1, allele2, kidallele) == expected
