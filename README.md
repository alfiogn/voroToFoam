# voroToFoam - Voronoi meshes for OpenFOAM #
A tool to generate Voronoi meshes for [OpenFOAM-10](https://github.com/OpenFOAM/OpenFOAM-10)
coupled with [voro++](https://math.lbl.gov/voro++/).

NOTE: it should be compatible also with older versions of OpenFOAM. Tested on [OpenFOAM-7](https://github.com/OpenFOAM/OpenFOAM-10).

## Setup and installation ##
First install the voro++ library from [this site](https://math.lbl.gov/voro++/download/) and follow the installation instructions.

Once the underlying OpenFOAM version is loaded, run the `Allwmake` script to compile the application.

NOTE: please change the library path if needed in `voroToFoam/Make/options`
```-I<VORO++_LIBRARY_PATH>/voro++```

## Usage ##
A cloud of points must be available in a file with the following structure:
```
<label> <scalar> <scalar> <scalar>
id0 x0 y0 z0
id1 x1 y1 z1
...
```

A typical syntax can then be:
```
voroToFoam -tol 1e-8 -cubeSide 0.5 points.dat
```

It is advisable to run `renumberMesh -overwrite` after.

NOTE: I have used this tool together with [tetgen](https://wias-berlin.de/software/tetgen/) or [GMSH](https://gmsh.info/).

## Tutorials

### Cube
Starting from a cloud of points in `points.dat`, the following grid is obtained
<img src="./tutorials/cube/screenshot.png"  width="300" height="300">



## TODO ##
* Avoid zero area faces: [reference](applications/voroToFoam/voronoiMesh.C?plain=1#L304)
* The mesh generation can be extended to generic shapes, like in [this site](https://math.lbl.gov/voro++/examples/)
* Parallel implementation

## References ##

* [Negrini, G., Parolini, N., & Verani, M. (2023). *The Rhie-Chow stabilized Box Method for the Stokes problem*.](https://arxiv.org/abs/2308.01059)
```
@misc{negrini2023rhiechow,
      title = {The {R}hie-{C}how stabilized {B}ox {M}ethod for the {S}tokes problem},
      author = {G. Negrini and N. Parolini and M. Verani},
      year = {2023},
      eprint = {2308.01059},
      archivePrefix = {arXiv},
      primaryClass = {math.NA}
}
```
* [Negrini, G., Parolini, N., & Verani, M. (2021). *A diffuse interface box method for elliptic problems*. Applied Mathematics Letters, 120, 107314.](https://www.sciencedirect.com/science/article/pii/S0893965921001646)
```
@article{NEGRINI2021107314,
    title = {A diffuse interface box method for elliptic problems},
    journal = {Applied Mathematics Letters},
    volume = {120},
    pages = {107314},
    year = {2021},
    issn = {0893-9659},
    doi = {https://doi.org/10.1016/j.aml.2021.107314},
    url = {https://www.sciencedirect.com/science/article/pii/S0893965921001646},
    author = {G. Negrini and N. Parolini and M. Verani},
    keywords = {Box method, Diffuse interface, Complex geometries}
}
```
