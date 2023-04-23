import pytest

from nrpytutorial.dirtree import NodeKind, Node
from nrpytutorial.dirtree import make_node

@pytest.mark.parametrize('testcase', ['astring'])
def test_making_unsupported_kind_rises_exception(testcase):
    with pytest.raises(TypeError):
        make_node('name', testcase)
