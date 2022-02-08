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

#include "fvcLaplacianPlus.H"
#include "laplacianScheme2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    
template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
laplacianScheme2Grad
(
    const GeometricField<GType, fvsPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme2<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    ).ref().fvcLaplacianGrad(gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
laplacianScheme2Grad
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme2<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    ).ref().fvcLaplacianGrad(gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
laplacianScheme2Grad
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    GeometricField<GType, fvsPatchField, surfaceMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fvc::laplacianScheme2Grad(Gamma, vf, name);
}

template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
laplacianScheme2Grad
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh>>& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> Laplacian
    (
        fvc::laplacianScheme2Grad(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
