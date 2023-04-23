import pytest

from nrpytutorial.dirtree import Dirtree
from pathlib import Path

class TestBssnDirtreeStructure:

    @pytest.fixture
    def dirtree():
        pass

    def test_root(dirtree):
        assert dirtree.get_node('root') == Path('bssn')

    def test_codes(dirtree):
        assert dirtree.get_node('codes') == Path('bssn').joinpath('ccodes')

    def test_output(dirtree):
        assert dirtree.get_node('output') == Path('bssn').joinpath('ccodes/output')
