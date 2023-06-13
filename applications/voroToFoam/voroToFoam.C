/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    voroToFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "voronoiMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate Voronoi tessellation from a set of points.\n"
        "The tessellation is then bounded by a cube (default L=1.0).\n\n"
        "The point file should be structured as:\n"
        "id0 x0 y0 z0\n"
        "id1 x1 y1 z1\n"
        "..."
    );

    argList::validArgs.append("fileName");
    argList::addOption("cubeSide", "scalar", "length of the cube side");
    argList::addOption("tol", "scalar", "tolerance value for geometric search");

    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    word nodesFile = args.argRead<word>(1);
    scalar L = args.optionLookupOrDefault("cubeSide", 1.0);
    scalar tol = args.optionLookupOrDefault("tol", small);

    voronoiMesh voroMesh(nodesFile, L, tol);

    Info<< nl << "Creating polyMesh" << endl;

    word defaultFacesName = "defaultFaces";
    word defaultFacesType = "patch";

    //Info<< voroMesh.IDs() << endl;
    //Info<< voroMesh.owner() << endl;
    //Info<< voroMesh.neighbour() << endl;

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        voroMesh.points(),
        voroMesh.faces(),
        //voroMesh.cells()
        voroMesh.owner(),
        voroMesh.neighbour()
    );

    List<polyPatch*> boundary(1);

    boundary[0] = new polyPatch
    (
        defaultFacesName,
        voroMesh.boundarySize(),
        voroMesh.start(),
        0,
        mesh.boundaryMesh(),
        defaultFacesType
    );

    mesh.addPatches(boundary);

    mesh.write();

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
