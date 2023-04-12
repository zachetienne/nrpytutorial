import pytest
from nrpytutorial import sample
def test_scenario_1(): assert sample() == 0
def test_scenario_2(): assert not sample() != 0
