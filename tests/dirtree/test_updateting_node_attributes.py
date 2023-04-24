import pytest

from nrpytutorial.dirtree import make_dir_node, make_file_node

@pytest.fixture
def dir_node():
    return make_dir_node('dir')

@pytest.fixture
def file_node():
    return make_file_node('file')

def test_adding_child_into_dir_node(dir_node, file_node):
    dir_node.add_child(file_node)
    assert dir_node.children == [file_node]

def test_adding_child_into_file_node_raises_exception(dir_node, file_node):
    with pytest.raises(ValueError):
        file_node.add_child(dir_node)

def test_updating_parent_with_file_node_raises_exception(dir_node, file_node):
    with pytest.raises(ValueError):
        dir_node.update_parent(file_node)

def test_updating_node_parent(dir_node, file_node):
    file_node.update_parent(dir_node)
    assert file_node.parent == dir_node
