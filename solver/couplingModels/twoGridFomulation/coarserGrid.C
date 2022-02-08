/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "coarserGrid.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::coarserGrid::coarserGrid
(
    const fvMesh& mesh,
    const vector coarseGridOrigin,
    const vector coarseGridVertex,
    const scalar coarseGridSize       
)
:       
    mesh_(mesh),
    
    coarseGridOrigin_(coarseGridOrigin),
    
    coarseGridVertex_(coarseGridVertex),
    
    coarseGridSize_(coarseGridSize),
    
    fineGridCellCenters_(mesh.cellCentres()),
        
    fineGridCellVolumes_(mesh.cellVolumes())
    
{
        
    Info << "Construct caose grid"<< endl;
    
    //iteration index
    label i,j,k;
    
    //coarse grid geometry
    coarseGridNumberInX_ = ceil(mag(coarseGridVertex_.x()-coarseGridOrigin_.x())/coarseGridSize);
    coarseGridNumberInY_ = ceil(mag(coarseGridVertex_.y()-coarseGridOrigin_.y())/coarseGridSize);
    coarseGridNumberInZ_ = ceil(mag(coarseGridVertex_.z()-coarseGridOrigin_.z())/coarseGridSize);
    
    scalar halfCoarseGridSize = coarseGridSize/2.0;
    
    
    //creat coarse grid
    label coarseGridNumber_ = coarseGridNumberInX_*coarseGridNumberInY_*coarseGridNumberInZ_;
    
    DynamicList<vector> coarseGridCenters(0);
    DynamicList<label> coarseGridIDList(0);
    
    vector coarseGridCenter(0, 0, 0);
    label coarseGridID = 0;
    
    for (k=1;k<=coarseGridNumberInZ_;k++)
    {
        for (j=1;j<=coarseGridNumberInY_;j++)
            {
                for (i=1;i<=coarseGridNumberInX_;i++)
                {
                    coarseGridCenter.x() = coarseGridOrigin_.x() + halfCoarseGridSize + coarseGridSize*(i-1);
                    coarseGridCenter.y() = coarseGridOrigin_.y() + halfCoarseGridSize + coarseGridSize*(j-1);
                    coarseGridCenter.z() = coarseGridOrigin_.z() + halfCoarseGridSize + coarseGridSize*(k-1);
                                            
                    coarseGridIDList.append(coarseGridID);
                    coarseGridCenters.append(coarseGridCenter);
                    
                    coarseGridID++;
                }
            }
    }
    
    coarseGridCellCenters_ = coarseGridCenters;
    coarseGridIDList_ = coarseGridIDList;
    
    //map fine grid to coarse grid
    labelListList coarseGridToFineGrid(coarseGridNumber_);
    scalar fineCellX_;
    scalar fineCellY_;
    scalar fineCellZ_;
    
    scalar xl_, xr_;
    scalar yl_, yr_;
    scalar zl_, zr_;
    
    forAll(coarseGridCellCenters_, idC)
    {
        DynamicList<label> fineGridIDToACoarseGrid(0);
        
        xl_ = coarseGridCellCenters_[idC].x() - halfCoarseGridSize;
        xr_ = coarseGridCellCenters_[idC].x() + halfCoarseGridSize;
        yl_ = coarseGridCellCenters_[idC].y() - halfCoarseGridSize;
        yr_ = coarseGridCellCenters_[idC].y() + halfCoarseGridSize;
        zl_ = coarseGridCellCenters_[idC].z() - halfCoarseGridSize;
        zr_ = coarseGridCellCenters_[idC].z() + halfCoarseGridSize;
        
        forAll(fineGridCellCenters_, idF)
        {
            fineCellX_ = fineGridCellCenters_[idF].x();
            fineCellY_ = fineGridCellCenters_[idF].y();
            fineCellZ_ = fineGridCellCenters_[idF].z();
            
            if ( (xl_ <= fineCellX_) && (fineCellX_ < xr_) && (yl_ <= fineCellY_) && (fineCellY_ < yr_) && (zl_ <= fineCellZ_) && (fineCellZ_ < zr_) )
            {
                fineGridIDToACoarseGrid.append(idF);
            }
        }
        
        coarseGridToFineGrid[idC] = fineGridIDToACoarseGrid;
    }
    
    coarseGridToFineGrid_ = coarseGridToFineGrid;
}
    
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//average scalar fields with weight fields, default weighting fields is uniform
Foam::tmp<volScalarField>
Foam::coarserGrid::coarserGrid::averagedField
(
    const volScalarField fineField, 
    volScalarField weightField, 
    volScalarField phaseFractionField
) const
{
    tmp<volScalarField> averagedFieldTmp
    (
        volScalarField::New
        (
            "averagedField",
            mesh_,
            fineField.dimensions()
        )
    );

    volScalarField& field = averagedFieldTmp.ref();
    field = fineField;
    scalarField& fieldCells = field.primitiveFieldRef();
    
    const scalarField& fine = fineField;
    const scalarField& weight = weightField;
    const scalarField& alpha = phaseFractionField;
    
    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sum = 0.0;
        scalar average = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();
        
        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sum = sum + fine[finIdList[idF]]*weight[finIdList[idF]]*alpha[finIdList[idF]];
                sumWeight = sumWeight + weight[finIdList[idF]]*alpha[finIdList[idF]];
            }
            
            average = sum/sumWeight;
            
            forAll(finIdList, idF)
            {
                fieldCells[finIdList[idF]] = average;
            }
        }
        
    }
    
    return averagedFieldTmp;
}
        
//average vector fields with weight fields, default weighting fields is uniform tmp<volVectorField> averagedField(volVectorField fineField, volScalarField weightField, volScalarField phaseFractionField) const;
tmp<volVectorField>
Foam::coarserGrid::coarserGrid::averagedField
(
    const volVectorField fineField, 
    volScalarField weightField, 
    volScalarField phaseFractionField
) const
{
    tmp<volVectorField> averagedFieldTmp
    (
        volVectorField::New
        (
            "averagedField",
            mesh_,
            fineField.dimensions()
        )
    );

    volVectorField& field = averagedFieldTmp.ref();
    field = fineField;
    vectorField& fieldCells = field.primitiveFieldRef();
    
    const vectorField& fine = fineField;
    const scalarField& weight = weightField;
    const scalarField& alpha = phaseFractionField;
    
    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sumx = 0.0;
        scalar sumy = 0.0;
        scalar sumz = 0.0;
        scalar averagex = 0.0;
        scalar averagey = 0.0;
        scalar averagez = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();

        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sumx = sumx + fine[finIdList[idF]].x()*weight[finIdList[idF]]*alpha[finIdList[idF]];
                sumy = sumy + fine[finIdList[idF]].y()*weight[finIdList[idF]]*alpha[finIdList[idF]];
                sumz = sumz + fine[finIdList[idF]].z()*weight[finIdList[idF]]*alpha[finIdList[idF]];
                sumWeight = sumWeight + weight[finIdList[idF]]*alpha[finIdList[idF]];
            }
            
            averagex = sumx/sumWeight;
            averagey = sumy/sumWeight;
            averagez = sumz/sumWeight;
            
            forAll(finIdList, idF)
            {
                fieldCells[finIdList[idF]].x() = averagex;
                fieldCells[finIdList[idF]].y() = averagey;
                fieldCells[finIdList[idF]].z() = averagez;
            }
        }
        
    }
    
    return averagedFieldTmp;
}  

//average scalar fields weighted by cell volume
tmp<volScalarField>
Foam::coarserGrid::coarserGrid::averagedField
(
    const volScalarField fineField, 
    volScalarField phaseFractionField
) const
{
    tmp<volScalarField> averagedFieldTmp
    (
        volScalarField::New
        (
            "averagedField",
            mesh_,
            fineField.dimensions()
        )
    );

    volScalarField& field = averagedFieldTmp.ref();
    field = fineField;
    scalarField& fieldCells = field.primitiveFieldRef();
    
    const scalarField& fine = fineField;
    const scalarField& weight = fineGridCellVolumes_;
    const scalarField& alpha = phaseFractionField;
    
    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sum = 0.0;
        scalar average = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();
        
        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sum = sum + fine[finIdList[idF]]*weight[finIdList[idF]]*alpha[finIdList[idF]];
                sumWeight = sumWeight + weight[finIdList[idF]]*alpha[finIdList[idF]];
            }
            
            average = sum/sumWeight;

            forAll(finIdList, idF)
            {
                fieldCells[finIdList[idF]] = average;
            }
        }
        
    }

    return averagedFieldTmp;
}

//average vector fields weighted by cell volume
tmp<volVectorField>
Foam::coarserGrid::coarserGrid::averagedField
(
    const volVectorField fineField, 
    volScalarField phaseFractionField
) const
{
    tmp<volVectorField> averagedFieldTmp
    (
        volVectorField::New
        (
            "averagedField",
            mesh_,
            fineField.dimensions()
        )
    );

    volVectorField& field = averagedFieldTmp.ref();
    field = fineField;
    vectorField& fieldCells = field.primitiveFieldRef();
    
    const vectorField& fine = fineField;
    const scalarField& weight = fineGridCellVolumes_;
    const scalarField& alpha = phaseFractionField;

    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sumx = 0.0;
        scalar sumy = 0.0;
        scalar sumz = 0.0;
        scalar averagex = 0.0;
        scalar averagey = 0.0;
        scalar averagez = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();

        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sumx = sumx + fine[finIdList[idF]].x()*weight[finIdList[idF]]*alpha[finIdList[idF]];
                sumy = sumy + fine[finIdList[idF]].y()*weight[finIdList[idF]]*alpha[finIdList[idF]];
                sumz = sumz + fine[finIdList[idF]].z()*weight[finIdList[idF]]*alpha[finIdList[idF]];
                sumWeight = sumWeight + weight[finIdList[idF]]*alpha[finIdList[idF]];
            }
            
            averagex = sumx/sumWeight;
            averagey = sumy/sumWeight;
            averagez = sumz/sumWeight;
            
            forAll(finIdList, idF)
            {
                fieldCells[finIdList[idF]].x() = averagex;
                fieldCells[finIdList[idF]].y() = averagey;
                fieldCells[finIdList[idF]].z() = averagez;
            }
        }
        
    }
    
    return averagedFieldTmp;
} 
        
//average scalar fields weighted by cell volume
tmp<volScalarField>
Foam::coarserGrid::coarserGrid::averagedField
(
    const volScalarField fineField
) const
{
    tmp<volScalarField> averagedFieldTmp
    (
        volScalarField::New
        (
            "averagedField",
            mesh_,
            fineField.dimensions()
        )
    );

    volScalarField& field = averagedFieldTmp.ref();
    field = fineField;
    scalarField& fieldCells = field.primitiveFieldRef();
    
    const scalarField& fine = fineField;
    const scalarField& weight = fineGridCellVolumes_;
    
    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sum = 0.0;
        scalar average = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();
        
        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sum = sum + fine[finIdList[idF]]*weight[finIdList[idF]];
                sumWeight = sumWeight + weight[finIdList[idF]];
            }
            
            average = sum/sumWeight;
            
            forAll(finIdList, idF)
            {
                fieldCells[finIdList[idF]] = average;
            }
        }
    }
    
    return averagedFieldTmp;
}

//average scalar fields weighted by cell volume
tmp<volScalarField>
Foam::coarserGrid::coarserGrid::averagedAlphaField
(
    const volScalarField fineVpField
) const
{
    tmp<volScalarField> averagedFieldTmp
    (
        volScalarField::New
        (
            "averagedField",
            mesh_,
            dimless
        )
    );

    volScalarField& field = averagedFieldTmp.ref();
    field = fineVpField/dimensionedScalar(dimVolume, 1);
    scalarField& fieldCells = field.primitiveFieldRef();
    
    const scalarField& fine = fineVpField;
    const scalarField& weight = fineGridCellVolumes_;

    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sum = 0.0;
        scalar average = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();
        
        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sum = sum + fine[finIdList[idF]];
                sumWeight = sumWeight + weight[finIdList[idF]];
            }
            
            average =1 - sum/sumWeight;
            
            forAll(finIdList, idF)
            {
                fieldCells[finIdList[idF]] = average;
            }
        }
    }

    return averagedFieldTmp;
}

    
Foam::tmp<scalarField>
Foam::coarserGrid::coarserGrid::averagedSource
(
    const scalarField& fineField, 
    volScalarField weightField
) const
{
    tmp<scalarField> averagedFieldTmp(new scalarField(fineField.size()));
    
    scalarField& field = averagedFieldTmp.ref();
    
    const scalarField& fine = fineField;
    const scalarField& weight = weightField;
    const scalarField& volWeight = fineGridCellVolumes_;
    
    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sum = 0.0;
        scalar average = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();
        
        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sum = sum + fine[finIdList[idF]]*weight[finIdList[idF]]*volWeight[finIdList[idF]];
                sumWeight = sumWeight + weight[finIdList[idF]]*volWeight[finIdList[idF]];
            }
            
            average = sum/sumWeight;
            
            forAll(finIdList, idF)
            {
                field[finIdList[idF]] = average;
            }
        }
        
    }
    
    return averagedFieldTmp;
}
        
//average vector fields with weight fields, default weighting fields is uniform tmp<volVectorField> averagedField(volVectorField fineField, volScalarField weightField, volScalarField phaseFractionField) const;
tmp<vectorField>
Foam::coarserGrid::coarserGrid::averagedSource
(
    const vectorField& fineField, 
    volScalarField weightField
) const
{
    tmp<vectorField> averagedFieldTmp(new vectorField(fineField.size()));
    
    vectorField& field = averagedFieldTmp.ref();
    
    const vectorField& fine = fineField;
    const scalarField& weight = weightField;
    const scalarField& volWeight = fineGridCellVolumes_;
    
    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sumx = 0.0;
        scalar sumy = 0.0;
        scalar sumz = 0.0;
        scalar averagex = 0.0;
        scalar averagey = 0.0;
        scalar averagez = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();

        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sumx = sumx + fine[finIdList[idF]].x()*weight[finIdList[idF]]*volWeight[finIdList[idF]];
                sumy = sumy + fine[finIdList[idF]].y()*weight[finIdList[idF]]*volWeight[finIdList[idF]];
                sumz = sumz + fine[finIdList[idF]].z()*weight[finIdList[idF]]*volWeight[finIdList[idF]];
                sumWeight = sumWeight + weight[finIdList[idF]]*volWeight[finIdList[idF]];
            }
            
            averagex = sumx/sumWeight;
            averagey = sumy/sumWeight;
            averagez = sumz/sumWeight;
            
            forAll(finIdList, idF)
            {
                field[finIdList[idF]].x() = averagex;
                field[finIdList[idF]].y() = averagey;
                field[finIdList[idF]].z() = averagez;
            }
        }
        
    }
    
    return averagedFieldTmp;
}  
        
Foam::tmp<scalarField>
Foam::coarserGrid::coarserGrid::averagedSource
(
    const scalarField& fineField
) const
{
    tmp<scalarField> averagedFieldTmp(new scalarField(fineField.size()));
    
    scalarField& field = averagedFieldTmp.ref();
    
    const scalarField& fine = fineField;
    const scalarField& volWeight = fineGridCellVolumes_;
    
    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sum = 0.0;
        scalar average = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();
        
        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sum = sum + fine[finIdList[idF]]*volWeight[finIdList[idF]];
                sumWeight = sumWeight + volWeight[finIdList[idF]];
            }
            
            average = sum/sumWeight;
            
            forAll(finIdList, idF)
            {
                field[finIdList[idF]] = average;
            }
        }
        
    }
    
    return averagedFieldTmp;
}
        
//average vector fields with weight fields, default weighting fields is uniform tmp<volVectorField> averagedField(volVectorField fineField, volScalarField weightField, volScalarField phaseFractionField) const;
tmp<vectorField>
Foam::coarserGrid::coarserGrid::averagedSource
(
    const vectorField& fineField
) const
{
    tmp<vectorField> averagedFieldTmp(new vectorField(fineField.size()));
    
    vectorField& field = averagedFieldTmp.ref();
    
    const vectorField& fine = fineField;
    const scalarField& volWeight = fineGridCellVolumes_;
    
    forAll(coarseGridToFineGrid_, idC)
    {
        scalar sumx = 0.0;
        scalar sumy = 0.0;
        scalar sumz = 0.0;
        scalar averagex = 0.0;
        scalar averagey = 0.0;
        scalar averagez = 0.0;
        scalar sumWeight = 0.0;
        labelList finIdList = coarseGridToFineGrid_[idC];
        label fineCellNumber = finIdList.size();

        if (fineCellNumber > 0)
        {
            forAll(finIdList, idF)
            {
                sumx = sumx + fine[finIdList[idF]].x()*volWeight[finIdList[idF]];
                sumy = sumy + fine[finIdList[idF]].y()*volWeight[finIdList[idF]];
                sumz = sumz + fine[finIdList[idF]].z()*volWeight[finIdList[idF]];
                sumWeight = sumWeight + volWeight[finIdList[idF]];
            }
            
            averagex = sumx/sumWeight;
            averagey = sumy/sumWeight;
            averagez = sumz/sumWeight;
            
            forAll(finIdList, idF)
            {
                field[finIdList[idF]].x() = averagex;
                field[finIdList[idF]].y() = averagey;
                field[finIdList[idF]].z() = averagez;
            }
        }
        
    }
    
    return averagedFieldTmp;
}  

// ************************************************************************* //
