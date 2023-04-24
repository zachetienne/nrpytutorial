import pytest

from nrpytutorial.dirtree import Dirtree, make_dir_node
from pathlib import Path

@pytest.fixture
def bssn_root_node():
    return make_dir_node('bssn')

@pytest.fixture
def bssn_codes_node(bssn_root_node):
    return make_dir_node('codes')

@pytest.fixture
def bssn_dirtree(bssn_root_node, bssn_codes_node):
    dirtree = Dirtree()
    dirtree.set_root(bssn_root_node)
    dirtree.add_node(bssn_codes_node, 'codes')
    return dirtree

def test_adding_node_without_root_raises_exception():
    node = make_dir_node('asdf')
    dirtree = Dirtree()
    with pytest.raises(AttributeError):
        dirtree.add_node(node, 'adsf')

@pytest.mark.skip(reason='work in progres')
@pytest.mark.parametrize('testcase, alias, expected', [
    ('bssn_dirtree', 'root', 'bssn_root_node'),
    ('bssn_dirtree', 'codes', 'bssn_codes_node')
    ])
def test_getting_node(testcase, alias, expected, request):
    testcase = request.getfixturevalue(testcase)
    actual = testcase.get_node(alias)
    expected = request.getfixturevalue(expected)
    assert actual == expected
