import numpy as np
from numpy.testing import assert_raises
from pysplit import make_trajectorygroup


def test_bogus_signature():
    signature = './nothing_to_see_here*'
    assert_raises(LookupError, make_trajectorygroup, signature)


def test_bogus_filenames():
    filenames = ['./not_a_trajectory0',
                 './not_a_trajectory1',
                 './not_a_trajectory2']
    assert_raises(IOError, make_trajectorygroup, filenames)


if __name__ == '__main__':
    np.testing.run_module_suite()
