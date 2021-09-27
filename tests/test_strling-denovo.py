import sys
sys.path.append("..")
from denovo import *
import pytest

@pytest.mark.parametrize("allele, wgl, minwig, expected", [
    (100, 0.1, 1, (90, 110.00000000000001)),
    (0, 0.1, 1, (-1, 1)),
    (10, 0.1, 5, (5, 15)),
])
def test_wiggle(allele, wgl, minwig, expected):
    assert wiggle(allele, wgl, minwig) == expected

# wiggle should raise an error when wgl is not between 0 and 1
@pytest.mark.parametrize("allele, wgl, minwig", [
    (100, 10, 1),
    (100, -10, 1),
])
def test_wiggle_error(allele, wgl, minwig):
    with pytest.raises(ValueError):
        wiggle(allele, wgl, minwig)

# to do
#@pytest.mark.parametrize("allele, wgl, minwig, expected", [
#    (10, 0.1, (1,500), (1,500) ),
#])
#def test_check_range(minwig, proportion, allele1, allele2, kidallele, expected):
#    assert check_range(minwig, proportion, allele1, allele2, kidallele) == expected

