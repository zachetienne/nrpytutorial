import pytest

from nrpytutorial.dirtree import Dirtree
from nrpytutorial.dirtree import DirtreeNodeDirector
from pathlib import Path


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
    return director.construct_file_node('header.h', should_backup=True)

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

def test_dirtree_backups_nodes_when_creating(dirtree, tmp_path, header_file):
    dirtree.create()
    dirtree.create()
    header_path = tmp_path.joinpath(header_file.path)
    header_backup = f"{header_path}_bak0"
    assert Path(header_backup).exists()

def test_dirtree_creates_subnode_inplace(dirtree, tmp_path):
    dirtree.make_dir('name', 'codes')
    assert tmp_path.joinpath(dirtree.get('codes').path).joinpath('name').exists()
