/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "voro++.hh"

#include "SortableList.H"
#include "voronoiMesh.H"

// * * * * * * * * * * * *   Private Member Functions  * * * * * * * * * * * //

void Foam::voronoiMesh::generateVoro() const
{
    Info<< "    - Set domain parameters" << nl;
    const double x_min = -side_, x_max = side_;
    const double y_min = -side_, y_max = side_;
    const double z_min = -side_, z_max = side_;
    const double cvol = (x_max - x_min)*(y_max - y_min)*(x_max - x_min);

    // Set up the number of blocks that the container is divided into
    const int n_x = 6, n_y = 6, n_z = 6;


    Info<< "    - Generate Voronoi mesh from seed" << nl;

    // Initialize Voronoi mesh
    voro::container con
    (
        x_min, x_max, y_min, y_max, z_min, z_max,
        n_x, n_y, n_z,
        false, false, false, 8
    );

    // Add particles from file
    con.import(fileName_.c_str());

    // Get cell number
    nCells_ = con.total_particles();
    nFaces_ = 0;

    cellVertices_.setSize(nCells_);
    cellFaces_.setSize(nCells_);
    cellFaceCell_.setSize(nCells_);

    Info<< "    - Loop Voronoi cells and build connectivity" << nl;

    // Loop on Voronoi cells
    voro::c_loop_all vl(con);
    int ijk, q; double *pp;
    voro::voronoicell_neighbor c;
    if (vl.start())
    {
        do if (con.compute_cell(c, vl))
        {
            ijk = vl.ijk;
            q = vl.q;
            pp = con.p[ijk] + con.ps*q;

            // Cell barycenter
            double cx = *pp;
            double cy = pp[1];
            double cz = pp[2];

            // Get cell ID
            int ID = con.id[ijk][q];

            // Get cell vertices
            std::vector<double> vertices;
            c.vertices(cx, cy, cz, vertices);

            // Get cell face cardinality
            std::vector<int> faceOrders;
            c.face_orders(faceOrders);

            // Get face-vertic3es connectivity
            std::vector<int> faceVertices;
            c.face_vertices(faceVertices);

            // Get cell neighbors IDs
            std::vector<int> faceNeighs;
            c.neighbors(faceNeighs);


            // Fill mesh connectivity structures
            // cell vertices
            int nCoords = vertices.size();
            cellVertices_[ID].setSize(nCoords/3, Zero);
            int nVerts = 0;
            int i = 0;
            while (i < nCoords)
            {
                cellVertices_[ID][nVerts][0] = vertices[i++];
                cellVertices_[ID][nVerts][1] = vertices[i++];
                cellVertices_[ID][nVerts++][2] = vertices[i++];
            }

            // cell faces
            int nFaces = faceOrders.size();
            nFaces_ += nFaces;
            cellFaces_[ID].setSize(nFaces);
            int count = 0;
            for (int i = 0; i < nFaces; i++)
            {
                cellFaces_[ID][i].setSize(faceOrders[i]);

                count++;
                for (int j = 0; j < faceOrders[i]; j++)
                {
                    cellFaces_[ID][i][j] = faceVertices[count++];
                }
            }

            // cell neighbours
            int nNeighs = faceNeighs.size();
            cellFaceCell_[ID].setSize(nNeighs, 0);
            for (int i = 0; i < nNeighs; i++)
            {
                cellFaceCell_[ID][i] = faceNeighs[i];
            }

        } while(vl.inc());
    }
}


void Foam::voronoiMesh::generateConn() const
{
    cellPoints_.setSize(nCells_);
    forAll(cellVertices_, ci)
    {
        cellPoints_[ci].setSize(cellVertices_[ci].size(), 0);
    }

    Info<< "    - Gather all polyhedra vertices" << nl;
    DynamicList<vector> points0(0);
    DynamicList<label> cells0(0);
    DynamicList<label> cellPoints0(0);

    forAll(cellVertices_, ci)
    {
        points0.append(cellVertices_[ci]);
        for (int i = 0; i < cellVertices_[ci].size(); i++)
        {
            cells0.append(ci);
            cellPoints0.append(i);
        }
    }

    Info<< "    - Sorting by x component" << nl;
    SortableList<scalar> xSortedPoints
    (
        vectorField(points0).component(0)
    );
    const labelList& xOrder = xSortedPoints.indices();

    Info<< "    - Find unique vertices" << nl;
    label nPoints = points0.size();
    points_.setSize(nPoints);
    pointsCell_.setSize(nPoints);

    label nPts = 0;
    label markDifferentX = 0;
    forAll(xSortedPoints, pi)
    {
        bool isUnique = true;

        label pointIdx = nPts;

        bool checked = false;

        for (int j = markDifferentX; j < nPts; j++)
        {
            if (mag(points0[xOrder[pi]][0] - points_[j][0]) < tol_)
            {
                if (!checked)
                {
                    markDifferentX = j;
                    checked = true;
                }

                if
                (
                    (mag(points0[xOrder[pi]][1] - points_[j][1]) < tol_)
                 && (mag(points0[xOrder[pi]][2] - points_[j][2]) < tol_)
                )
                {
                    isUnique = false;
                    pointIdx = j;
                }
            }

            if (!isUnique)
            {
                break;
            }
        }

        cellPoints_
        [
            cells0[xOrder[pi]]
        ][
            cellPoints0[xOrder[pi]]
        ] = pointIdx;

        if (isUnique)
        {
            pointsCell_[nPts].first() = cells0[xOrder[pi]];
            pointsCell_[nPts].second() = cellPoints0[xOrder[pi]];
            points_[nPts++] = points0[xOrder[pi]];
        }
    }

    pointsCell_.setSize(nPts);
    points_.setSize(nPts);

    Info<<"       Original points = " << points0.size() << nl;
    Info<<"       Unique points   = " << nPts << nl << endl;
}


void Foam::voronoiMesh::generateMesh() const
{
    cells_.setSize(nCells_, cell(0));

    Info<< "    - Find unique faces" << nl;
    List<labelList> faces(nFaces_);
    label nBoundaryFaces = 0;
    DynamicList<label> own(0);
    DynamicList<label> nei(0);

    label nFaces = 0;
    forAll(cellFaces_, ci)
    {
        const labelListList& facesI = cellFaces_[ci];
        cells_[ci].setSize(facesI.size(), -1);

        forAll(facesI, fj)
        {
            label faceIndex = -1;

            const labelList& faceJ = facesI[fj];

            // check face is unique
            bool isUnique = true;

            if (cellFaceCell_[ci][fj] < 0)
            {
                nBoundaryFaces++;
            }
            else
            {
                const cell& cellJ = cells_[cellFaceCell_[ci][fj]];

                if (cellJ.size())
                {
                    // We passed this cell already
                    isUnique = false;

                    for (int k = 0; k < cellJ.size(); k++)
                    {
                        if (own[cellJ[k]] == ci)
                        {
                            faceIndex = cellJ[k];
                            break;
                        }
                    }
                }
            }

            if (isUnique)
            {
                // register face
                nei.append(ci);
                own.append(cellFaceCell_[ci][fj]);

                faces[nFaces].setSize(faceJ.size(), -1);

                forAll(faceJ, pk)
                {
                    faces[nFaces][pk] = cellPoints_[ci][faceJ[pk]];
                }

                //if (face(faces[nFaces]).mag(points_) < small)
                //{
                //    Info<< "found face with zero area"<<nl;
                //}

                faceIndex = nFaces;

                nFaces++;
            }

            cells_[ci][fj] = faceIndex;
        }
    }

    faces.setSize(nFaces);

    nFaces_ = nFaces;


    Info<< "    - Check boundary faces and finalize" << nl;
    faces_.setSize(faces.size());
    owner_.setSize(own.size(), -1);
    neighbour_.setSize(nei.size(), -1);

    faceList boundaryFaces(faces.size());
    labelList boundaryOwn(own.size(), -1);

    start_ = faces.size() - nBoundaryFaces;
    bdSize_ = nBoundaryFaces;

    label internalFaces = 0;
    nBoundaryFaces = 0;

    forAll(faces, fi)
    {
        if (own[fi] < 0)
        {
            boundaryOwn[nBoundaryFaces] = nei[fi];
            boundaryFaces[nBoundaryFaces].setSize(faces[fi].size());

            for (int i=0; i<faces[fi].size();i++)
            {
                boundaryFaces[nBoundaryFaces][i] =
                    faces[fi][faces[fi].size() - i - 1];
            }
            nBoundaryFaces++;
        }
        else
        {
            owner_[internalFaces] = own[fi];
            neighbour_[internalFaces] = nei[fi];
            faces_[internalFaces++] = face(faces[fi]);
        }
    }

    owner_.setSize(internalFaces);
    neighbour_.setSize(internalFaces);
    faces_.setSize(internalFaces);

    boundaryFaces.setSize(nBoundaryFaces);
    boundaryOwn.setSize(nBoundaryFaces);
    faces_.append(boundaryFaces);
    owner_.append(boundaryOwn);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::voronoiMesh::voronoiMesh
(
    const word& fileName,
    const scalar side,
    const scalar tol
)
:
    fileName_(fileName),
    side_(side),
    tol_(tol),
    nCells_(-1),
    nFaces_(-1),
    cellVertices_(),
    cellFaces_(),
    cellFaceCell_(),
    seedPoints_(),
    points_(),
    cells_(),
    faces_(),
    owner_(),
    neighbour_(),
    bdSize_(-1),
    start_(-1),
    pointsCell_(),
    cellPoints_()
{
    Info<< nl << "Generating voronoi mesh from " << fileName << nl;
    generateVoro();

    Info<< nl << "Generating connectivity" << nl;
    generateConn();

    Info<< nl << "Generating mesh components" << nl;
    generateMesh();
}


// ************************************************************************* //
