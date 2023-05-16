import pytest

from nrpytutorial.dirtree import _make_dir_node, _make_file_node

@pytest.fixture
def parent(tmp_path):
    return _make_dir_node(str(tmp_path))

@pytest.fixture
def codes(parent):
    return _make_dir_node('codes', parent=parent)

@pytest.fixture
def lib(parent):
    return _make_file_node('lib', parent=parent)

@pytest.mark.parametrize('testcase', ['codes', 'lib'])
def test_node_path_conversion(testcase, tmp_path, request):
    actual = request.getfixturevalue(testcase)
    expected = tmp_path.joinpath(testcase)
    assert actual.to_path() == expected

@pytest.mark.parametrize('testcase', ['codes', 'lib'])
def test_node_string_conversion(testcase, tmp_path, request):
    actual = request.getfixturevalue(testcase)
    expected = tmp_path.joinpath(testcase)
    assert actual.to_string() == str(expected)
