import pytest

from nrpytutorial.dirtree import Dirtree
from pathlib import Path

@pytest.fixture
def bssn_nodes():
    return {'root': Path('bssn'),
            'codes': Path('bssn').joinpath('codes'),
            'build': Path('bssn').joinpath('build'),
            'output': Path('bssn').joinpath('build/output')
            }


# class TestBssnDirtreeStructure:
#
#     @pytest.fixture
#     def dirtree():
#         pass
#
#     def test_root(dirtree):
#         assert dirtree.get_node('root') == Path('bssn')
#
#     def test_codes(dirtree):
#         assert dirtree.get_node('codes') == Path('bssn').joinpath('ccodes')
#
#     def test_output(dirtree):
#         assert dirtree.get_node('output') == Path(
#             'bssn').joinpath('ccodes/output')
