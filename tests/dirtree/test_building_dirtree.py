import pytest

from nrpytutorial.dirtree import Dirtree, make_dir_node
from pathlib import Path

@pytest.fixture
def bssn_root_node():
    return make_dir_node('bssn')

@pytest.fixture
def bssn_dirtree(bssn_root_node):
    dirtree = Dirtree()
    dirtree.set_root(bssn_root_node)
    return dirtree

def test_adding_node_without_root_raises_exception():
    node = make_dir_node('asdf')
    dirtree = Dirtree()
    with pytest.raises(AttributeError):
        dirtree.add_node(node, 'adsf')

def test_getting_node_using_alias(bssn_dirtree, bssn_root_node):
    actual = bssn_dirtree.get_node('root')
    assert actual == bssn_root_node
