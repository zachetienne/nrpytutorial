import enum


@enum.unique
class CoordSystemKind(enum.Enum):
    SPHERICAL = 'Spherical'
    SINHSPHERICAL = "SinhSpherical"
    SINHSPHERICALV2 = 'SinhSphericalv2'
    CYLINDRICAL = "Cylindrical"
    SINHCYLINDRICAL = "SinhCylindrical"
    SYMTP = "SymTP"
    SINHSYMTP = "SinhSymTP"

    def as_string(self):
        return self.value
