import enum
from dataclasses import dataclass
from functools import partial


@enum.unique
class PulseProfile(enum.Enum):
    GAUSSIAN = 'Gaussian_pulse'
    GAUSSIAN_V2 = 'Gaussian_pulsev2'
    TANH = 'Tanh_pulse'

    def as_string(self):
        return self.value


@dataclass
class PulseParameters:
    profile: PulseProfile
    amplitude: float
    center: float
    width: float

    def as_tuple(self):
        return (self.profile.as_string(),
                self.amplitude,
                self.center,
                self.width,
                )


def build_pulse(profile: PulseProfile,
                amplitude: float,
                center: float,
                width: float):
    return PulseParameters(profile, amplitude, center, width)


build_gaussian_pulse = partial(build_pulse, PulseProfile.GAUSSIAN)
