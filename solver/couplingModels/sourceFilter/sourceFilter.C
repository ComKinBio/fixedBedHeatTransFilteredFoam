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

#include "sourceFilter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::sourceFilter::sourceFilter
(
    word cloudName,
    const fvMesh& mesh,
    Time& diffusionRunTime,
    const fvMesh& diffusionmesh,
    simpleControl& simple
)
:       
    cloudName_(cloudName),
    mesh_(mesh),
    diffusionRunTime_(diffusionRunTime),
    diffusionMesh_(diffusionmesh),
    simple_(simple),
    particleProperties_
    (
        IOobject
        (
            cloudName_ + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    constantProperties_
    (
        particleProperties_.subDict("constantProperties")
    ),
    diffusionBandWidth_(constantProperties_, "diffusionBandWidth", 0.024),
    diffusionBandWidthForMassCoupling_(constantProperties_, "diffusionBandWidthMass", 0.024),
    diffusionBandWidthForMomentumCoupling_(constantProperties_, "diffusionBandWidth", 0.024),
    diffusionBandWidthForHeatCoupling_(constantProperties_, "diffusionBandWidthMomentum", 0.024),
    diffusionBandWidthForRadiaCoupling_(constantProperties_, "diffusionBandWidthRadia", 0.024),
    diffusionSteps_(constantProperties_, "diffusionSteps", 6),
    implicitFvm_(constantProperties_, "useImplicitLaplacian", true),
    smoothDirection_
    (
        constantProperties_.lookupOrDefault
        (
            "smoothDirection",
            tensor(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0)
        )
    ),
    DT("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_),
    startTime(diffusionRunTime_.startTime()),
    startTimeIndex(diffusionRunTime_.startTimeIndex()),
    diffusionTime_(0),
    diffusionTimeForMass_(0),
    diffusionTimeForMomentum_(0),
    diffusionTimeForHeat_(0),
    diffusionTimeForRadia_(0),
    diffusionDeltaT_(0),
    diffusionDeltaTForMass_(0),
    diffusionDeltaTForMomentum_(0),
    diffusionDeltaTForHeat_(0),
    diffusionDeltaTForRadia_(0),
    phaseDiffusionFlag(false),
    massDiffusionFlag(false),
    momentumDiffusionFlag(false),
    heatDiffusionFlag(false),
    radiaDiffusionFlag(false)
{
    
    //- initialization information output
    if (!constantProperties_.found("diffusionBandWidth"))
    {
        Info<< "Diffusion band width for phase set to 0.024 by default" << endl;
    }
    else
    {
        Info<< "*** Diffusion band width for phase set to " << diffusionBandWidth_.value() << endl;
    }
    
    if (!constantProperties_.found("diffusionBandWidthMass"))
    {
        Info<< "Diffusion band width For Mass Coupling set to 0.024 by default" << endl;
    }
    else
    {
        Info<< "*** Diffusion band width For Mass Coupling set to " << diffusionBandWidthForMassCoupling_.value() << endl;
    }
    
    if (!constantProperties_.found("diffusionBandWidthMomentum"))
    {
        Info<< "Diffusion band width For Momentum Coupling set to 0.024 by default" << endl;
    }
    else
    {
        Info<< "*** Diffusion band width For Momentum Coupling set to " << diffusionBandWidthForMomentumCoupling_.value() << endl;
    }

    if (!constantProperties_.found("diffusionBandWidthHeat"))
    {
        Info<< "Diffusion band width For Heat Coupling set to 0.024 by default" << endl;
    }
    else
    {
        Info<< "*** Diffusion band width For Heat Coupling set to " << diffusionBandWidthForHeatCoupling_.value() << endl;
    }
    
    if (!constantProperties_.found("diffusionBandWidthRadia"))
    {
        Info<< "Diffusion band width For Radia Coupling set to 0.024 by default" << endl;
    }
    else
    {
        Info<< "*** Diffusion band width For Radia Coupling set to " << diffusionBandWidthForRadiaCoupling_.value() << endl;
    }
    
    if (!constantProperties_.found("diffusionSteps"))
    {
        Info<< "*** Diffusion steps set to 6 by default"<< endl;
    }
    else
    {
        Info<< "*** Diffusion steps set to " << diffusionSteps_.value() << endl;
    }
    
    
    // determine the time and time step in diffusion procedure
    diffusionTime_ = pow(diffusionBandWidth_.value(), 2)/4;
    diffusionDeltaT_ = diffusionTime_/diffusionSteps_.value();

    // determine the time and time step for mass in diffusion procedure
    diffusionTimeForMass_ = pow(diffusionBandWidthForMassCoupling_.value(), 2)/4;
    diffusionDeltaTForMass_ = diffusionTimeForMass_/diffusionSteps_.value();

    // determine the time and time step for heat in diffusion procedure
    diffusionTimeForHeat_ = pow(diffusionBandWidthForHeatCoupling_.value(), 2)/4;
    diffusionDeltaTForHeat_ = diffusionTimeForHeat_/diffusionSteps_.value();

    // determine the time and time step for momentum in diffusion procedure
    diffusionTimeForMomentum_ = pow(diffusionBandWidthForMomentumCoupling_.value(), 2)/4;
    diffusionDeltaTForMomentum_ = diffusionTimeForMomentum_/diffusionSteps_.value();

    // determine the time and time step for Radia in diffusion procedure
    diffusionTimeForRadia_ = pow(diffusionBandWidthForRadiaCoupling_.value(), 2)/4;
    diffusionDeltaTForRadia_ = diffusionTimeForRadia_/diffusionSteps_.value();
    
    //- set diffusion flags for different source terms
    if (diffusionBandWidth_.value() > 0)
    {
        phaseDiffusionFlag = true;
    }
    else
    {
        phaseDiffusionFlag = false;
    }
    
    if (diffusionBandWidthForMassCoupling_.value() > 0)
    {
        massDiffusionFlag = true;
    }
    else
    {
        massDiffusionFlag = false;
    }
    
    if (diffusionBandWidthForMomentumCoupling_.value() > 0)
    {
        momentumDiffusionFlag = true;
    }
    else
    {
        momentumDiffusionFlag = false;
    }
    
    if (diffusionBandWidthForHeatCoupling_.value() > 0)
    {
        heatDiffusionFlag = true;
    }
    else
    {
        heatDiffusionFlag = false;
    }
    
    if (diffusionBandWidthForRadiaCoupling_.value() > 0)
    {
        radiaDiffusionFlag = true;
    }
    else
    {
        radiaDiffusionFlag = false;
    }
    
}
    
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline Foam::scalar
Foam::sourceFilter::diffusionBandWidth() const
{
    return diffusionBandWidth_.value();
}


inline Foam::scalar
Foam::sourceFilter::diffusionBandWidthForMassCoupling() const
{
    return diffusionBandWidthForMassCoupling_.value();
}


inline Foam::scalar
Foam::sourceFilter::diffusionBandWidthForMomentumCoupling() const
{
    return diffusionBandWidthForMomentumCoupling_.value();
}


inline Foam::scalar
Foam::sourceFilter::diffusionBandWidthForHeatCoupling() const
{
    return diffusionBandWidthForHeatCoupling_.value();
}


inline Foam::scalar
Foam::sourceFilter::diffusionBandWidthForRadiaCoupling() const
{
    return diffusionBandWidthForRadiaCoupling_.value();
}


inline Foam::label
Foam::sourceFilter::diffusionSteps() const
{
    return diffusionSteps_.value();
}


bool Foam::sourceFilter::getPhaseDiffusionFlag() const
{
    return phaseDiffusionFlag;
}


bool Foam::sourceFilter::getMassDiffusionFlag() const
{
    return massDiffusionFlag;
}


bool Foam::sourceFilter::getMomentumDiffusionFlag() const
{
    return momentumDiffusionFlag;
}


bool Foam::sourceFilter::getHeatDiffusionFlag() const
{
    return heatDiffusionFlag;
}


bool Foam::sourceFilter::getRadiaDiffusionFlag() const
{
    return radiaDiffusionFlag;
}

Foam::scalar Foam::sourceFilter::diffusionEndTime(word type)
{
    if (type == "phase")
    {
        return diffusionTime_;
    }
    else if (type == "mass")
    {
        return diffusionTimeForMass_;
    }
    else if (type == "momentum")
    {
        return diffusionTimeForMomentum_;
    }
    else if (type == "heat")
    {
        return diffusionTimeForHeat_;
    }
    else if (type == "radia")
    {
        return diffusionTimeForRadia_;
    }
    else
    {
        FatalErrorInFunction
            << "source type is not defined" << type
            << abort(FatalError);
    }
    
    return 0;
}
    
        
Foam::scalar Foam::sourceFilter::diffusionDeltaT(word type)
{
    if (type == "phase")
    {
        return diffusionDeltaT_;
    }
    else if (type == "mass")
    {
        return diffusionDeltaTForMass_;
    }
    else if (type == "momentum")
    {
        return diffusionDeltaTForMomentum_;
    }
    else if (type == "heat")
    {
        return diffusionDeltaTForHeat_;
    }
    else if (type == "radia")
    {
        return diffusionDeltaTForRadia_;
    }
    else
    {
        FatalErrorInFunction
            << "source type is not defined" << type
            << abort(FatalError);
    }
    
     return 0;
}
        
void Foam::sourceFilter::updateDiffusionTime()
{
    diffusionTime_ = pow(diffusionBandWidth_.value(), 2)/4;
    diffusionDeltaT_ = diffusionTime_/diffusionSteps_.value();

    // determine the time and time step for mass in diffusion procedure
    diffusionTimeForMass_ = pow(diffusionBandWidthForMassCoupling_.value(), 2)/4;
    diffusionDeltaTForMass_ = diffusionTimeForMass_/diffusionSteps_.value();

    // determine the time and time step for heat in diffusion procedure
    diffusionTimeForHeat_ = pow(diffusionBandWidthForHeatCoupling_.value(), 2)/4;
    diffusionDeltaTForHeat_ = diffusionTimeForHeat_/diffusionSteps_.value();

    // determine the time and time step for momentum in diffusion procedure
    diffusionTimeForMomentum_ = pow(diffusionBandWidthForMomentumCoupling_.value(), 2)/4;
    diffusionDeltaTForMomentum_ = diffusionTimeForMomentum_/diffusionSteps_.value();

    // determine the time and time step for Radia in diffusion procedure
    diffusionTimeForRadia_ = pow(diffusionBandWidthForRadiaCoupling_.value(), 2)/4;
    diffusionDeltaTForRadia_ = diffusionTimeForRadia_/diffusionSteps_.value();
}

void Foam::sourceFilter::diffusion
(
    volScalarField& s,
    word type
)
{
    
    scalar endTime = diffusionEndTime(type);
    scalar deltaT = diffusionDeltaT(type);
    
    diffusionRunTime_.setEndTime(endTime);
    diffusionRunTime_.setDeltaT(deltaT);
    
    volScalarField diffWorkField
    (
        IOobject
        (
            "tempDiffScalar",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedScalar
        (
            "zero",
            dimless,
            scalar(0.0)
        ),
        zeroGradientFvPatchScalarField::typeName
        
    );
    
//     diffWorkField.internalField() = s.internalField();
    scalarField& diffWorkFieldInterFeildRef = diffWorkField.ref();
    
    scalarField& sInterFeildRef = s.ref();

    diffWorkFieldInterFeildRef = sInterFeildRef;
    
    if (implicitFvm_.value())
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
            }
        }

    }
    else
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT, diffWorkField));
            }
        }
    }
    
//     s.internalField() = diffWorkField.internalField();
    sInterFeildRef = diffWorkFieldInterFeildRef;
    
    diffusionRunTime_.setTime(startTime, startTimeIndex);

    return;
}
        

Foam::tmp<Foam::volScalarField::Internal> Foam::sourceFilter::diffusion
(
    const volScalarField::Internal& s,
    word type
)
{    
    
    scalar endTime = diffusionEndTime(type);
    scalar deltaT = diffusionDeltaT(type);
    
    diffusionRunTime_.setEndTime(endTime);
    diffusionRunTime_.setDeltaT(deltaT);
    
    tmp<volScalarField::Internal> tS
    (
        volScalarField::Internal::New
        (
            "tS",
            mesh_,
            dimensionedScalar(s.dimensions(), 0.0)
        )
    );

    scalarField& S = tS.ref();
    
    S = s;
    
    volScalarField diffWorkField
    (
        IOobject
        (
            "tempDiffScalar",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedScalar
        (
            "zero",
            s.dimensions(),
            scalar(0.0)
        ),
        zeroGradientFvPatchScalarField::typeName
        
    );

    
    scalarField& diffWorkFieldInterFeildRef = diffWorkField.ref();
    
    diffWorkFieldInterFeildRef = S;
    
    if (implicitFvm_.value())
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
            }
        }

    }
    else
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT, diffWorkField));
            }
        }
    }
    
    S = diffWorkField.internalField();
    
    diffusionRunTime_.setTime(startTime, startTimeIndex);

    return tS;
}

Foam::tmp<Foam::volVectorField::Internal> Foam::sourceFilter::diffusion
(
    const volVectorField::Internal& s,
    word type
)
{
    
    scalar endTime = diffusionEndTime(type);
    scalar deltaT = diffusionDeltaT(type);
    
    diffusionRunTime_.setEndTime(endTime);
    diffusionRunTime_.setDeltaT(deltaT);
    
    tmp<volVectorField::Internal> tS
    (
        volVectorField::Internal::New
        (
            "tS",
            mesh_,
            dimensionedVector(s.dimensions(), vector::zero)
        )
    );

    vectorField& S = tS.ref();
    
    S = s;
    
    volVectorField diffWorkField
    (
        IOobject
        (
            "tempDiffVector",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedVector
        (
            "zero",
            s.dimensions(),
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
        
    );

    
    vectorField& diffWorkFieldInterFeildRef = diffWorkField.ref();
    
    diffWorkFieldInterFeildRef = S;
    
    if (implicitFvm_.value())
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
            }
        }

    }
    else
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT, diffWorkField));
            }
        }
    }
    
    S = diffWorkField.internalField();
    
    diffusionRunTime_.setTime(startTime, startTimeIndex);

    return tS;
}







// ************************************************************************* //
