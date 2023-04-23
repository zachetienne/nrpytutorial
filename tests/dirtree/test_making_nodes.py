import pytest

from nrpytutorial.dirtree import NodeKind, Node
from nrpytutorial.dirtree import make_node, make_file_node
from collections import namedtuple

@pytest.mark.parametrize('testcase', ['astring'])
def test_making_unsupported_kind_node_rises_exception(testcase):
    with pytest.raises(TypeError):
        make_node('name', testcase)

@pytest.mark.parametrize('factory, args', [
    (make_node, ('name', NodeKind.FILE)),
    (make_file_node, ('name',))
    ])
def test_making_file_node(factory, args):
    actual = factory(*args)
    assert actual.kind is NodeKind.FILE
