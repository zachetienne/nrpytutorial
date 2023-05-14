import pytest

from nrpytutorial.dirtree import make_dir_node, make_file_node
from nrpytutorial.dirtree import Dirtree
from pathlib import Path

@pytest.fixture
def bssn_root_dir():
    return make_dir_node('bssn')

@pytest.fixture
def bssn_codes_dir():
    return make_dir_node('codes')

@pytest.fixture
def bssn_main_file():
    return make_file_node('main.c')

@pytest.fixture
def bssn_dirtree(bssn_root_dir, bssn_codes_dir, bssn_main_file):
    dirtree = Dirtree()
    dirtree.update_root(bssn_root_dir)
    dirtree.add_node(bssn_codes_dir, 'root')
    dirtree.add_node(bssn_main_file, 'codes', key='main')
    return dirtree

def test_adding_node_without_root_raises_exception():
    node = make_dir_node('asdf')
    dirtree = Dirtree()
    with pytest.raises(AttributeError):
        dirtree.add_node(node, 'adsf')

@pytest.mark.parametrize('testcase, key, expected', [
    ('bssn_dirtree', 'root', 'bssn_root_dir'),
    ('bssn_dirtree', 'codes', 'bssn_codes_dir'),
    ('bssn_dirtree', 'main', 'bssn_main_file')
    ])
def test_getting_node(testcase, key, expected, request):
    testcase = request.getfixturevalue(testcase)
    actual = testcase.get_node(key)
    expected = request.getfixturevalue(expected)
    assert actual == expected

def test_getting_node_path(bssn_dirtree):
    actual = bssn_dirtree.get_node_path('main')
    assert actual == Path('bssn').joinpath('codes/main.c')
