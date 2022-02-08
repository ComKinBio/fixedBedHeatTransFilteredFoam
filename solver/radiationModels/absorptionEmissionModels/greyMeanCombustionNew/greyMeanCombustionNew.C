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

\*---------------------------------------------------------------------------*/

#include "greyMeanCombustionNew.H"
#include "combustionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(greyMeanCombustionNew, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        greyMeanCombustionNew,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMeanCombustionNew::
greyMeanCombustionNew
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    greyMean(dict, mesh, typeName),
    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMeanCombustionNew::
~greyMeanCombustionNew()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMeanCombustionNew::ECont
(
    const label bandI
) const
{
    tmp<volScalarField> E = greyMean::ECont(bandI);

    const word& name = combustionModel::combustionPropertiesName;
    E.ref() += EhrrCoeff_*mesh_.lookupObject<combustionModel>(name).Qdot();

    return E;
}

void Foam::radiationModels::absorptionEmissionModels::greyMeanCombustionNew::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aj
) const
{
   a = this->a();
   
   forAll(aj, bandI)
   {
       aj[bandI] =  a;
   }

}


// ************************************************************************* //
