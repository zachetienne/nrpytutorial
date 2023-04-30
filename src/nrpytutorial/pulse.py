import enum
from functools import partial
from collections import namedtuple


@enum.unique
class PulseProfile(enum.Enum):
    GAUSSIAN = 'Gaussian_pulse'
    GAUSSIAN_V2 = 'Gaussian_pulsev2'
    TANH = 'Tanh_pulse'

    def as_string(self):
        return self.value

PulseParameters = namedtuple('PulseParameters', 'profile amplitude center width')

def build_pulse(profile: PulseProfile,
                amplitude: float,
                center: float,
                width: float):
    return PulseParameters(profile, amplitude, center, width)


build_gaussian_pulse = partial(build_pulse, PulseProfile.GAUSSIAN)
