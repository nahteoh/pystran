# pystran: Python for Structural Analysis

A simple structural analysis tool in Python.
The intent is for this package to be used in a classroom setting,
and the students are expected to fill in missing functionality.
Lot of the source code is included, such as for plotting,
but the basic mechanics formulations are not implemented on purpose.

![Alt pystran capabilities in graphic abstract](docs/splash.png "pystran")

## News

- 01/22/2025: Implement initial functionality. 

[Past news](#past-news)

## Limitations

- Two-dimensional and three-dimensional structures made up of truss (axial)
  members and beams (even in combination) can be handled.
- Only elastic models can be solved.
- Member loading is not considered. All member loading needs to be converted to nodal forces.
- Internal hinges can be modelled with linked joints.
- Offsets are currently not implemented.

## Requirements

- NumPy
- SciPy
- Matplotlib

These requirements can be easily satisfied by running the examples in Spyder 6 (IDE).

## Running

This package is not distributed through the official Python channels.
It needs to be downloaded from GitHub as a zip file, and expanded in some convenient location. 

The easiest way to run a pystran example is to download and install Spyder 6. It is a complete IDE for Python, including a very capable debugger. Just open an example and click the run button.

It is also possible to run using a plain Python.
The user then needs to install the requirements, and in the
pystran folder rather an example for instance as
```
py examples/linked_cantilevers_prescribed.py
```

## <a name="past-news"></a>Past news
