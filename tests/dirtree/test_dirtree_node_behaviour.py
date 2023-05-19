import pytest

from nrpytutorial.dirtree import DirtreeNodeDirector
from pathlib import Path

@pytest.fixture
def director():
    return DirtreeNodeDirector()

@pytest.fixture
def root_node(director, tmp_path):
    return director.construct_dir_node(str(tmp_path))

@pytest.fixture
def file_node(director, root_node):
    return director.construct_file_node('name', parent=root_node)

@pytest.fixture
def dir_node(director, root_node):
    return director.construct_dir_node('name', parent=root_node)

@pytest.fixture
def combined_node(file_node, dir_node):
    file_node.parent = dir_node
    file_node.data.should_backup = True
    return file_node

@pytest.mark.parametrize('testcase, expected', [
    ['file_node', 'name'],
    ['dir_node', 'name'],
    ['combined_node', 'name/name']
    ])
def test_node_path_representation(testcase, expected, request, tmp_path):
    actual = request.getfixturevalue(testcase)
    assert actual.path == tmp_path.joinpath(expected)

@pytest.mark.parametrize('testcase, expected', [
    ['file_node', 'name'],
    ['dir_node', 'name'],
    ['combined_node', 'name/name']
    ])
def test_node_string_representation(testcase, expected, request, tmp_path):
    actual = request.getfixturevalue(testcase)
    expected = tmp_path.joinpath(expected)
    assert str(actual) == str(expected)

@pytest.mark.parametrize('testcase, expected', [
    ['file_node', 'name'],
    ['dir_node', 'name'],
    ['combined_node', 'name/name']
    ])
def test_node_is_created(tmp_path, testcase, expected, request):
    actual = request.getfixturevalue(testcase)
    actual.create()
    assert tmp_path.joinpath(expected).exists()

def test_node_is_backuped(tmp_path, combined_node):
    combined_node.create()
    combined_node.backup()
    assert tmp_path.joinpath('name/name_v0').exists()
