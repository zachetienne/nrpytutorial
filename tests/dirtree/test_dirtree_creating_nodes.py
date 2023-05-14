import pytest

from nrpytutorial.dirtree import Dirtree

@pytest.fixture
def bssn_root():
    return 'build/bssn'

@pytest.fixture
def bssn_dirtree(tmp_path, bssn_root):
    root = tmp_path.joinpath(bssn_root)
    dirtree = Dirtree(root)
    return dirtree

@pytest.fixture
def create_bssn_once(bssn_dirtree):
    bssn_dirtree.create()

@pytest.mark.parametrize('expected', ['bssn_root'])
def test_dirtree_creates_nodes(tmp_path, create_bssn_once, expected, request):
    expected = request.getfixturevalue(expected)
    expected = tmp_path.joinpath(expected)
    assert expected.exists()
