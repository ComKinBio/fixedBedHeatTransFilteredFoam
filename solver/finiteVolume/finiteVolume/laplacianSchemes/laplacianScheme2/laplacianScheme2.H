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

Class
    Foam::fv::laplacianScheme2

Description
    Abstract base class for laplacian schemes.

SourceFiles
    laplacianScheme2.C

\*---------------------------------------------------------------------------*/

#ifndef laplacianScheme2_H
#define laplacianScheme2_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "linear.H"
#include "correctedSnGrad.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class fvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class laplacianScheme2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class laplacianScheme2
:
    public tmp<laplacianScheme2<Type, GType>>::refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;
        tmp<surfaceInterpolationScheme<GType>> tinterpGammaScheme_;
        tmp<snGradScheme<Type>> tsnGradScheme_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construction
        laplacianScheme2(const laplacianScheme2&);

        //- Disallow default bitwise assignment
        void operator=(const laplacianScheme2&) = delete;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            laplacianScheme2,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        laplacianScheme2(const fvMesh& mesh)
        :
            mesh_(mesh),
            tinterpGammaScheme_(new linear<GType>(mesh)),
            tsnGradScheme_(new correctedSnGrad<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        laplacianScheme2(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh),
            tinterpGammaScheme_(nullptr),
            tsnGradScheme_(nullptr)
        {
            tinterpGammaScheme_ = tmp<surfaceInterpolationScheme<GType>>
            (
                surfaceInterpolationScheme<GType>::New(mesh, is)
            );

            tsnGradScheme_ = tmp<snGradScheme<Type>>
            (
                snGradScheme<Type>::New(mesh, is)
            );
        }

        //- Construct from mesh, interpolation and snGradScheme schemes
        laplacianScheme2
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            mesh_(mesh),
            tinterpGammaScheme_(igs),
            tsnGradScheme_(sngs)
        {}


    // Selectors

        //- Return a pointer to a new laplacianScheme2 created on freestore
        static tmp<laplacianScheme2<Type, GType>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~laplacianScheme2();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
        
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fvcLaplacianGrad
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
        
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fvcLaplacianGrad
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) =0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvLaplacianTypeScheme2(SS, GType, Type)                             \
    typedef Foam::fv::SS<Foam::Type, Foam::GType> SS##Type##GType;             \
    defineNamedTemplateTypeNameAndDebug(SS##Type##GType, 0);                   \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            typedef SS<Type, GType> SS##Type##GType;                           \
                                                                               \
            laplacianScheme2<Type, GType>::                                     \
                addIstreamConstructorToTable<SS<Type, GType>>                  \
                add##SS##Type##GType##IstreamConstructorToTable_;              \
        }                                                                      \
    }


#define makeFvLaplacianScheme2(SS)                                              \
                                                                               \
makeFvLaplacianTypeScheme2(SS, scalar, scalar)                                  \
makeFvLaplacianTypeScheme2(SS, symmTensor, scalar)                              \
makeFvLaplacianTypeScheme2(SS, tensor, scalar)                                  \
makeFvLaplacianTypeScheme2(SS, scalar, vector)                                  \
makeFvLaplacianTypeScheme2(SS, symmTensor, vector)                              \
makeFvLaplacianTypeScheme2(SS, tensor, vector)                                  \
makeFvLaplacianTypeScheme2(SS, scalar, sphericalTensor)                         \
makeFvLaplacianTypeScheme2(SS, symmTensor, sphericalTensor)                     \
makeFvLaplacianTypeScheme2(SS, tensor, sphericalTensor)                         \
makeFvLaplacianTypeScheme2(SS, scalar, symmTensor)                              \
makeFvLaplacianTypeScheme2(SS, symmTensor, symmTensor)                          \
makeFvLaplacianTypeScheme2(SS, tensor, symmTensor)                              \
makeFvLaplacianTypeScheme2(SS, scalar, tensor)                                  \
makeFvLaplacianTypeScheme2(SS, symmTensor, tensor)                              \
makeFvLaplacianTypeScheme2(SS, tensor, tensor)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "laplacianScheme2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
