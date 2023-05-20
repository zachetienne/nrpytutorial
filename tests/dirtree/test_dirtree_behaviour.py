import pytest

from nrpytutorial.dirtree import Dirtree
from nrpytutorial.dirtree import DirtreeNodeDirector


@pytest.fixture
def director():
    return DirtreeNodeDirector()

@pytest.fixture
def bssn_node(director: DirtreeNodeDirector):
    return director.construct_dir_node('bssn')

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
    dirtree.add_node(lib_file, 'codes')
    dirtree.add_node(headers, 'codes')
    dirtree.add_node(header_file, 'headers')
    dirtree.create()
    return dirtree

def test_dirtree(dirtree):
    assert False
