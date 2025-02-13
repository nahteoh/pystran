# pystran: Python for Structural Analysis

A simple structural analysis tool in Python, for structures consisting of truss and beam members, both in two dimensions and in three dimensions.


![Alt pystran capabilities in graphic abstract](docs/splash.png)

## News

- 02/12/2025: Make it possible to use general identifiers.
- 02/05/2025: Add general springs to ground.
- 01/30/2025: Add tutorials.
- 01/22/2025: Implement initial functionality. 

[Past news](#past-news)

## Features & Limitations

- The package analyzes two-dimensional and three-dimensional structures made up
  of truss (axial) members and beams (even in combination), with added masses
  and springs at joints.
- Linear statics and dynamics (free vibration) solvers are included.
- The Bernoulli-Euler model is implemented, so no shear deformation is taken into account.
- Only elastic models can be solved.
- Only straight members are treated.
- Only doubly symmetric cross sections can be handled in three dimensions. Hence
  there is no coupling between the bending actions in the two orthogonal planes.
- Warping of the cross sections is not modelled, hence only free torsion effects are included.
- Member loading is not considered. All member loading needs to be converted to nodal forces.
- Internal hinges can be modelled with linked joints. No member end releases are implemented.
- Degrees of freedom are only along the cartesian axes. Skew supports are not
  included (except with a penalty method based on springs)
- Offsets are currently not implemented.

## Requirements

- NumPy
- SciPy
- Matplotlib

These requirements can be easily satisfied by running the examples in [Spyder 6](https://www.spyder-ide.org/download/) (IDE).

## Running

This package is not distributed through the official Python channels.
It needs to be downloaded from GitHub as a zip file, and expanded in some convenient location. 

The __pystran folder__ can be located by looking for this README.md file.

The easiest way to run a pystran example is to download and install [Spyder 6](https://www.spyder-ide.org/download/). It is a complete IDE for Python, 
including a very capable debugger. Just open an example and click the run button.

It is also possible to run simulations using a plain Python.
The user then needs to install the requirements, and in the
pystran folder run an example for instance as
```
py examples/linked_cantilevers_prescribed.py
```

## Tutorials

Step-by-step tutorials are available in the [`tutorials`](./tutorials) folder. 
Run tutorials for example as
```
py tutorials/three_bars_tut.py
```

## Examples

There are many examples in the [`examples`](./examples) folder. They may not be heavily documented,
but they do show many recipes for solving truss and beam structures.

## Documentation

The functions in the package are documented. The documentation in nicely formatted form can be obtained by running `pdoc`. In the `pystran` folder, run
```
pdoc pystran/ --math -o html
```
which will create an `./html` folder. double click the `index.html` file in that folder.

## <a name="past-news"></a>Past news
