import pytest

from nrpytutorial.dirtree import DirtreeNodeDirector
from pathlib import Path

@pytest.fixture
def director():
    return DirtreeNodeDirector()

@pytest.fixture
def file_node(director):
    return director.construct_file_node('name')

@pytest.fixture
def dir_node(director):
    return director.construct_dir_node('name')

@pytest.fixture
def combined_node(file_node, dir_node):
    file_node.parent = dir_node
    return file_node

@pytest.mark.parametrize('testcase, expected', [
    ['file_node', 'name'],
    ['dir_node', 'name'],
    ['combined_node', 'name/name']
    ])
def test_node_path_representation(testcase, expected, request):
    actual = request.getfixturevalue(testcase)
    assert actual.to_path() == Path(expected)

@pytest.mark.parametrize('testcase, expected', [
    ['file_node', 'name'],
    ['dir_node', 'name'],
    ['combined_node', 'name/name']
    ])
def test_node_string_representation(testcase, expected, request):
    actual = request.getfixturevalue(testcase)
    assert actual.to_string() == expected
