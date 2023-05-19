import pytest

from nrpytutorial.dirtree import DirtreeNodeDirector
from pathlib import Path

@pytest.fixture
def director():
    return DirtreeNodeDirector()

@pytest.mark.parametrize('testcase, expected',[
    (('create_file_node', ['name']), 'name'),
    (('create_dir_node', ['name']), 'name'),
    ])
def test_node_path_representation(director, testcase, expected):
    mthd, args = testcase
    mthd = getattr(director, mthd)
    actual = mthd(*args)
    assert actual.to_path() == Path(expected)

@pytest.mark.parametrize('testcase, expected',[
    (('create_file_node', ['name']), 'name'),
    (('create_dir_node', ['name']), 'name'),
    ])
def test_node_path_representation(director, testcase, expected):
    mthd, args = testcase
    mthd = getattr(director, mthd)
    actual = mthd(*args)
    assert actual.to_string() == expected
