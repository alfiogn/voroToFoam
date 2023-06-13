/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
\*---------------------------------------------------------------------------*/

#ifndef voronoiMesh_H
#define voronoiMesh_H

#include "typeInfo.H"
#include "labelList.H"
#include "Field.H"
#include "vectorList.H"
#include "faceList.H"
#include "cellList.H"
#include "labelPair.H"
#include "IFstream.H"
#include "IOstream.H"
#include "DynamicList.H"
#include "IStringStream.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class voronoiMesh;

/*---------------------------------------------------------------------------*\
                        Class voronoiMesh Declaration
\*---------------------------------------------------------------------------*/

class voronoiMesh
{
    // Private data

        //- Name of the seed file
        word fileName_;

        //- Cube side
        scalar side_;

        //- Tolerance to compare vertices
        scalar tol_;

        //- Cells number
        mutable label nCells_;

        //- Number of faces
        mutable label nFaces_;

        //- Cells labeling
        mutable labelList IDs_;

        //- Vertices of a cell
        mutable List<List<vector> > cellVertices_;

        //- Faces of a cell
        mutable List<labelListList> cellFaces_;

        //- Neighbours of a face
        mutable labelListList cellFaceCell_;

        //- Points of the mesh
        mutable vectorField seedPoints_;

        //- Points of the mesh
        mutable vectorField points_;

        //- Cells of the mesh
        mutable cellList cells_;

        //- Faces of the mesh
        mutable faceList faces_;

        //- Owner cells
        mutable labelList owner_;

        //- Neighbour cells
        mutable labelList neighbour_;

        //- Size of boundary faces
        mutable label bdSize_;

        //- First label of boundary faces
        mutable label start_;

        //- Other connectivity data
        mutable List<labelPair> pointsCell_;
        mutable List<List<label>> cellPoints_;


    // Private Member Functions

        //- Read file
        void generateVoro() const;

        //- Create connectivity
        void generateConn() const;

        //- Read file
        void generateMesh() const;


public:

    // Constructors

        //- Constructor
        voronoiMesh
        (
            const word&,
            const scalar side = 1.0,
            const scalar tol = small
        );

        //- Copy constructur disallowed
        voronoiMesh(const voronoiMesh&) = delete;


    // Destructor
    virtual ~voronoiMesh(){};


    // Public Member Functions

        //- Get points
        pointField points() const
        {
            return points_;
        }

        //- Get cells
        cellList cells() const
        {
            return cells_;
        }

        //- Get faces
        faceList faces() const
        {
            return faces_;
        }

        //- Get owner
        labelList owner() const
        {
            return owner_;
        }

        //- Get neighbour
        labelList neighbour() const
        {
            return neighbour_;
        }

        //- Get boundary size
        const label& boundarySize() const
        {
            return bdSize_;
        }

        //- Get boundary start face
        const label& start() const
        {
            return start_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
