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

#include "laplacianScheme2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define the constructor function hash tables

#define makeLaplacianGTypeScheme2(Type, GType)                                  \
    typedef laplacianScheme2<Type, GType> laplacianScheme2##Type##GType;         \
    defineTemplateRunTimeSelectionTable(laplacianScheme2##Type##GType, Istream);

#define makeLaplacianScheme2(Type)                                              \
    makeLaplacianGTypeScheme2(Type, scalar);                                    \
    makeLaplacianGTypeScheme2(Type, symmTensor);                                \
    makeLaplacianGTypeScheme2(Type, tensor);

makeLaplacianScheme2(scalar);
makeLaplacianScheme2(vector);
makeLaplacianScheme2(sphericalTensor);
makeLaplacianScheme2(symmTensor);
makeLaplacianScheme2(tensor);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
