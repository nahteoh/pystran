"""
pystran - Python package for structural analysis with trusses and beams 

(C) 2025, Petr Krysl

This package is intended for educational purposes only. Professional users may
find it too bare bones.

The approach is based on classical weighted residual formulation (Galerkin). The
formulations are derived in the "Finite element modeling with shells and beams"
[book](http://hogwarts.ucsd.edu/~pkrysl/femstructures-book/).


## Features and limitations

- Two-dimensional and three-dimensional structures made up of truss (axial)
  members and beams (even in combination) can be handled.
- The Bernoulli-Euler model is implemented, so no shear deformation is taken into account.
- Only elastic models can be solved.
- Only straight members are treated.
- Only doubly symmetric cross sections can be handled in three dimensions. Hence
  there is no coupling between the bending actions in the two orthogonal planes.
- Warping of the cross sections is not modelled, hence only free torsion effects are included.
- Member loading is not considered. All member loading needs to be converted to nodal forces.
- Internal hinges can be modelled with linked joints. No member end releases are implemented.
- Degrees of freedom are only along the cartesian axes. Skew supports are not included.
- Offsets are currently not implemented.

"""
