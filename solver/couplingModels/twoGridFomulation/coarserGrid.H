/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::coarserGrid

Description

SourceFiles
    coarserGrid.H
    coarserGrid.C

\*---------------------------------------------------------------------------*/

#ifndef coarserGrid_H
#define coarserGrid_H

#include "fvCFD.H"
#include "labelList.H"
#include "vectorList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coarserGrid Declaration
\*---------------------------------------------------------------------------*/

class coarserGrid
// public IOdictionary
{  
    //*******************************access and information read from settings*********************************//
//     //Thermo class
//     const psiReactionThermo& thermo_;
    
    //fine Mesh class
    const fvMesh& mesh_;
    
    //Two point to decide coarse grid geometry, the origin point
    const vector coarseGridOrigin_;
    
    //The second point to determine coarse grid geometry
    const vector coarseGridVertex_;
    
    //coarse grid size
    const scalar coarseGridSize_;
    
    //*******************************variables in grid calculation*********************************//

        //- fine grid center list 
        const vectorField& fineGridCellCenters_;
        
        //- fine grid volume list 
        const scalarField& fineGridCellVolumes_;
 
        //- coarse grid number in x direction
        label coarseGridNumberInX_;
        
        //- coarse grid number in y direction
        label coarseGridNumberInY_;
        
        //- coarse grid number in z direction
        label coarseGridNumberInZ_;
        
        //- coarse grid center list 
        vectorList coarseGridCellCenters_;
        
        //- coarse grid ID list
        labelList coarseGridIDList_;
        
        //- coarse Grid id To Fine Grid id
        labelListList coarseGridToFineGrid_;
        
        

    public:
  
    // Constructors

        //- Construct from fine mesh and input settings
        coarserGrid
        (
            const fvMesh& mesh,
            const vector coarseGridOrigin_,
            const vector coarseGridVertex_,
            const scalar coarseGridSize_       
        );


    //- Destructor
    ~coarserGrid()
    {}

        
    //average scalar fields with weight fields, default weighting fields is uniform
    tmp<volScalarField> averagedField(const volScalarField fineField, volScalarField weightField, volScalarField phaseFractionField) const;
        
    //average vector fields with weight fields, default weighting fields is uniform
    tmp<volVectorField> averagedField(const volVectorField fineField, volScalarField weightField, volScalarField phaseFractionField) const;
        
    //average scalar fields weighted by cell volume
    tmp<volScalarField> averagedField(const volScalarField fineField, volScalarField phaseFractionField) const;
        
    //average vector fields weighted by cell volume
    tmp<volVectorField> averagedField(const volVectorField fineField, volScalarField phaseFractionField) const;
        
    //average scalar fields weighted by cell volume
    tmp<volScalarField> averagedField(const volScalarField fineField) const;
    
    //average scalar fields weighted by cell volume
    tmp<volScalarField> averagedAlphaField(const volScalarField fineVpField) const;

    
            
    //average scalar fields with weight fields, default weighting fields is uniform
    tmp<scalarField> averagedSource(const scalarField& fineField, volScalarField weightField) const;
        
    //average vector fields with weight fields, default weighting fields is uniform
    tmp<vectorField> averagedSource(const vectorField& fineField, volScalarField weightField) const;
        
    //average scalar fields weighted by cell volume
    tmp<scalarField> averagedSource(const scalarField& fineField) const;
        
    //average vector fields weighted by cell volume
    tmp<vectorField> averagedSource(const vectorField& fineField) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
