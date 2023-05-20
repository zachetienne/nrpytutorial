import pytest

from nrpytutorial.dirtree import Dirtree
from nrpytutorial.dirtree import DirtreeNodeDirector


@pytest.fixture
def director():
    return DirtreeNodeDirector()

@pytest.fixture
def bssn_node(tmp_path, director: DirtreeNodeDirector):
    test_root = director.construct_dir_node(str(tmp_path))
    return director.construct_dir_node('bssn', parent=test_root)

@pytest.fixture
def codes_node(director: DirtreeNodeDirector):
    return director.construct_dir_node('codes')

@pytest.fixture
def lib_file(director: DirtreeNodeDirector):
    return director.construct_file_node('lib.c')

@pytest.fixture
def headers(director: DirtreeNodeDirector):
    return director.construct_dir_node('headers')

@pytest.fixture
def header_file(director: DirtreeNodeDirector):
    return director.construct_file_node('header.h')

@pytest.fixture
def dirtree(bssn_node, codes_node, lib_file, headers, header_file):
    dirtree = Dirtree(bssn_node)
    dirtree.add_codes(codes_node)
    dirtree.add(lib_file, 'codes')
    dirtree.add(headers, 'codes')
    dirtree.add(header_file, 'headers')
    return dirtree

@pytest.mark.parametrize('expected', ['bssn_node', 'codes_node', 'lib_file', 'headers', 'header_file'])
def test_dirtree_created_nodes(dirtree, expected, tmp_path, request):
    dirtree.create()
    expected = request.getfixturevalue(expected)
    assert tmp_path.joinpath(expected.path).exists()
