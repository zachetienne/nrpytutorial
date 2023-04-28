import pytest

from nrpytutorial.ccode_rendering import Cparameter, Cprototype
from nrpytutorial.ccode_rendering import CparameterFactory, CprototypeFactory


@pytest.fixture
def cparameter_factory():
    return CparameterFactory()


@pytest.fixture
def cprototype_factory(cparameter_factory):
    return CprototypeFactory(cparameter_factory)


@pytest.fixture
def ricci_eval_proto(cprototype_factory):
    return cprototype_factory.create_ricci_eval()


def test_cfunction_proto_rendering(ricci_eval_proto):
    actual = ricci_eval_proto.as_string()
    expected = "Ricci_eval(&rfmstruct, &params, RK_INPUT_GFS, auxevol_gfs);"
    assert actual == expected
