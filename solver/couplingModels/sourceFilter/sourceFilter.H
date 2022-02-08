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
    Foam::sourceFilter

Description

    filter method for particle source terms

SourceFiles
    sourceFilter.H
    sourceFilter.C

\*---------------------------------------------------------------------------*/

#ifndef sourceFilter_H
#define sourceFilter_H

#include "fvCFD.H"
#include "simpleControl.H"
#include "error.H"
#include "demandDrivenEntry.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coarserGrid Declaration
\*---------------------------------------------------------------------------*/

class sourceFilter
{  
    //*******************************access and information read from settings*********************************//
    //- cloud name
    word cloudName_;
    
    //- Gas Mesh class
    const fvMesh& mesh_;

    //- Diffusion time class
    Time& diffusionRunTime_;
    
    //- Diffusion control class
    const fvMesh& diffusionMesh_;
    
    //- Diffusion Mesh class
    simpleControl& simple_;
    
    //- Dictionary of particle properties
    IOdictionary particleProperties_;
    
    //- ConstantProperties dictionary
    const dictionary& constantProperties_;
    
    //- Phase diffusion band width
    demandDrivenEntry<scalar> diffusionBandWidth_;
    
    //- Mass diffusion band width
    demandDrivenEntry<scalar> diffusionBandWidthForMassCoupling_;
    
    //- Momentum diffusion band width
    demandDrivenEntry<scalar> diffusionBandWidthForMomentumCoupling_;
    
    //- Heat diffusion band width
    demandDrivenEntry<scalar> diffusionBandWidthForHeatCoupling_;
    
    //- Radiation diffusion band width
    demandDrivenEntry<scalar> diffusionBandWidthForRadiaCoupling_;
    
    //- Diffusion number of time steps
    demandDrivenEntry<label> diffusionSteps_;
    
    //- Scheme True if use implicit laplacian, False if use explicit laplacian
    demandDrivenEntry<bool> implicitFvm_;
    
    //- Diffusion tensor direction
    tensor smoothDirection_;
    
    //- Diffusion tensor
    dimensionedTensor DT;
    
    //- Diffusion start time
    dimensionedScalar startTime;
    label startTimeIndex;
   
    //- Diffusion time and time step for different source terms
    scalar diffusionTime_;
    scalar diffusionTimeForMass_;
    scalar diffusionTimeForMomentum_;
    scalar diffusionTimeForHeat_;
    scalar diffusionTimeForRadia_;
    
    scalar diffusionDeltaT_;
    scalar diffusionDeltaTForMass_;
    scalar diffusionDeltaTForMomentum_;
    scalar diffusionDeltaTForHeat_;
    scalar diffusionDeltaTForRadia_;
    

    //- flag whether the diffusion is applied to the source terms or phase
    bool phaseDiffusionFlag;
    bool massDiffusionFlag;
    bool momentumDiffusionFlag;
    bool heatDiffusionFlag;
    bool radiaDiffusionFlag;   
    
        
    
    //*******************************variables in grid calculation*********************************//      
        

    public:
  
    // Constructors

        //- Construct from fine mesh and input settings
        sourceFilter
        (
            word cloudName,
            const fvMesh& mesh,
            Time& diffusionRunTime,
            const fvMesh& diffusionmesh,
            simpleControl& simple
        );


    //- Destructor
    ~sourceFilter()
    {}

    // Member Functions
    
    
        //- Return const access to the Phase diffusion band width
        inline scalar diffusionBandWidth() const;
        
        //- Return const access to the Mass diffusion band width
        inline scalar diffusionBandWidthForMassCoupling() const;
        
        //- Return const access to the Momentum diffusion band width
        inline scalar diffusionBandWidthForMomentumCoupling() const;
        
        //- Return const access to the Heat diffusion band width
        inline scalar diffusionBandWidthForHeatCoupling() const;
        
        //- Return const access to the Radiation diffusion band width
        inline scalar diffusionBandWidthForRadiaCoupling() const;
        
        //- Return const access to the Diffusion number of time steps
        inline label diffusionSteps() const;
        
        
        //- Return const access to the flags of the source terms or phase
        bool getPhaseDiffusionFlag() const;
        
        bool getMassDiffusionFlag() const;
        
        bool getMomentumDiffusionFlag() const;
        
        bool getHeatDiffusionFlag() const;
        
        bool getRadiaDiffusionFlag() const;
        
        //- Get end time and dt according to source type
        //- the input shoule be "phase" "mass" "momentum" "heat" "radia"
        scalar diffusionEndTime(word type);
        
        scalar diffusionDeltaT(word type);
        
        //- update diffusion time from new diffusionBandWidth
        void updateDiffusionTime();
        
        //- Return the diffused field obtained from the given volScalarField s
        void diffusion(volScalarField& s, word type);
        
        //- Return the diffused field obtained from the given volScalarField s
        tmp<volScalarField::Internal> diffusion(const volScalarField::Internal& s, word type);

        //- Return the diffused field obtained from the given volVectorField s
        tmp<volVectorField::Internal> diffusion(const volVectorField::Internal& s, word type);
        
    

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
