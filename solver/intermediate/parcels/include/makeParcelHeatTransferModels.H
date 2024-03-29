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

#ifndef makeParcelHeatTransferModels_H
#define makeParcelHeatTransferModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NoHeatTransfer.H"
#include "RanzMarshall.H"
#include "KTA.H"
#include "fixedValue.H"
#include "fixedValue2.H"
#include "Wakao.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelHeatTransferModels(CloudType)                                \
                                                                               \
    makeHeatTransferModel(CloudType);                                          \
                                                                               \
    makeHeatTransferModelType(NoHeatTransfer, CloudType);                      \
    makeHeatTransferModelType(RanzMarshall, CloudType);                        \
    makeHeatTransferModelType(KTA, CloudType);                                 \
    makeHeatTransferModelType(fixedValue, CloudType);                          \
    makeHeatTransferModelType(fixedValue2, CloudType);                         \
    makeHeatTransferModelType(Wakao, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
