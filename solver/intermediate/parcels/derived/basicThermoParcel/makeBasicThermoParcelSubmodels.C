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
    modified by Tian
\*---------------------------------------------------------------------------*/

#include "basicThermoCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeParcelInjectionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeThermoParcelSurfaceFilmModels.H" // thermo variant

//Collision
#include "makeParcelCollisionModels.H"

//Collision Heat
#include "makeParcelCollisionHeatModels.H"

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(basicThermoCloud);

// Kinematic sub-models
makeThermoParcelForces(basicThermoCloud);
makeParcelDispersionModels(basicThermoCloud);
makeParcelInjectionModels(basicThermoCloud);
makeParcelPatchInteractionModels(basicThermoCloud);
makeParcelStochasticCollisionModels(basicThermoCloud);
makeParcelSurfaceFilmModels(basicThermoCloud);

//Collision sub-models
makeParcelCollisionModels(basicThermoCloud); 
makeParcelCollisionHeatModels(basicThermoCloud); 

// Thermo sub-models
makeParcelHeatTransferModels(basicThermoCloud);


// ************************************************************************* //
