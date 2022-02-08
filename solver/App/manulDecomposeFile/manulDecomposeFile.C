/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    manulDecomposeFile

Description
    This utility also read gridCoarseDict to generate decomposition file
    One processer size is required by the time of coarseGridsize

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "labelList.H"
#include "vectorList.H"
#include "DynamicList.H"
#include "IOobjectList.H"
#include "labelIOField.H"
#include "labelFieldIOField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //read input settings
    IOdictionary couplingDict
    (
        IOobject
        (
            "couplingDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    const dictionary TFMDict = couplingDict.subDict("TFM");

    vector gridOrigin;
    TFMDict.lookup("gridOrigin") >> gridOrigin;
    
    vector gridVertex;
    TFMDict.lookup("gridVertex") >> gridVertex;
    
    scalar gridSize = readScalar(TFMDict.lookup("gridSize"));
    
    //how many times that a proc is larger than the coarseGridSize
    label procSizeX = readLabel(TFMDict.lookup("procSizeX"));
    label procSizeY = readLabel(TFMDict.lookup("procSizeY"));
    label procSizeZ = readLabel(TFMDict.lookup("procSizeZ"));
    
    
    //*******************************variables initializtion*********************************//
    
    //Two point to decide Proc geometry, the origin point
    const vector procOrigin_ = gridOrigin;
    
    //The second point to determine Proc geometry
    const vector procVertex_ = gridVertex;
    
    //Proc size
    const scalar procSize_ = gridSize;
    
    //- fine grid center list 
    const vectorField& fineGridCellCenters_ = mesh.cellCentres();
        
    //- Proc number in x direction
    label procNumberInX_;
        
    //- Proc number in y direction
    label procNumberInY_;
        
    //- Proc number in z direction
    label procNumberInZ_;
        
    //- Proc center list 
    vectorList procCellCenters_;
        
    //- Proc ID list
    labelList procIDList_;

        
    //*******************************variables in grid calculation*********************************//
    
    //iteration index
    label i,j,k;
    
    //Proc geometry
    procNumberInX_ = ceil(mag(procVertex_.x()-procOrigin_.x())/(procSize_*procSizeX));
    procNumberInY_ = ceil(mag(procVertex_.y()-procOrigin_.y())/(procSize_*procSizeY));
    procNumberInZ_ = ceil(mag(procVertex_.z()-procOrigin_.z())/(procSize_*procSizeZ));
    
Info<<"procNuberInX: "<<procNumberInX_<<endl;
Info<<"procNuberInY: "<<procNumberInY_<<endl;
Info<<"procNuberInZ: "<<procNumberInZ_<<endl;
    
    scalar halfprocSizeX = procSize_*procSizeX/2.0;
    scalar halfprocSizeY = procSize_*procSizeY/2.0;
    scalar halfprocSizeZ = procSize_*procSizeZ/2.0;
    
    DynamicList<vector> procCenters(0);
    DynamicList<label> procIDList(0);
    
    vector procCenter(0, 0, 0);
    label procID = 0;
    
    for (k=1;k<=procNumberInZ_;k++)
    {
        for (j=1;j<=procNumberInY_;j++)
            {
                for (i=1;i<=procNumberInX_;i++)
                {
                    procCenter.x() = procOrigin_.x() + halfprocSizeX + procSize_*procSizeX*(i-1);
                    procCenter.y() = procOrigin_.y() + halfprocSizeY + procSize_*procSizeY*(j-1);
                    procCenter.z() = procOrigin_.z() + halfprocSizeZ + procSize_*procSizeZ*(k-1);
                                            
                    procIDList.append(procID);
                    procCenters.append(procCenter);
                    
                    procID++;
                }
            }
    }
    
    procCellCenters_ = procCenters;
    procIDList_ = procIDList;
    
    //map fine grid to Proc
    scalar fineCellX_;
    scalar fineCellY_;
    scalar fineCellZ_;
    
    scalar xl_, xr_;
    scalar yl_, yr_;
    scalar zl_, zr_;
    
    DynamicList<label> procIDListForFineGrid(0);
    
        
    forAll(fineGridCellCenters_, idF)
    {

        fineCellX_ = fineGridCellCenters_[idF].x();
        fineCellY_ = fineGridCellCenters_[idF].y();
        fineCellZ_ = fineGridCellCenters_[idF].z();
        
        label find = 0;
        label nofind = 0;
        
        forAll(procCellCenters_, idC)
        {
            xl_ = procCellCenters_[idC].x() - halfprocSizeX;
            xr_ = procCellCenters_[idC].x() + halfprocSizeX;
            yl_ = procCellCenters_[idC].y() - halfprocSizeY;
            yr_ = procCellCenters_[idC].y() + halfprocSizeY;
            zl_ = procCellCenters_[idC].z() - halfprocSizeZ;
            zr_ = procCellCenters_[idC].z() + halfprocSizeZ;
        
            if ( (xl_ <= fineCellX_) && (fineCellX_ < xr_) && (yl_ <= fineCellY_) && (fineCellY_ < yr_) && (zl_ <= fineCellZ_) && (fineCellZ_ < zr_)  && (find == 0) )
            {
                find++;
                procIDListForFineGrid.append(idC);
            }
        }
        
    }
    labelList procIds(fineGridCellCenters_.size());
            
    procIds = procIDListForFineGrid;
    
    labelIOList cellDecomposition
    (
        IOobject
        (
                "cellDecomposition",
                runTime.constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
        procIds
    );
    cellDecomposition.write();

    Info<< nl << "Wrote decomposition to "
        << cellDecomposition.objectPath()
        << " for use in manual decomposition." << endl;
                    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
