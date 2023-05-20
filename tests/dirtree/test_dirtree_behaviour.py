import pytest

from nrpytutorial.dirtree import Dirtree
from nrpytutorial.dirtree import DirtreeNodeDirector


@pytest.fixutre
def director():
    return DirtreeNodeDirector()

@pytest.fixture
def bssn(director: DirtreeNodeDirector):
    return director.construct_dir_node()
