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
    Foam::radiationModels::P1P

Description
    Particulate Effects in the P-1 Model

SourceFiles
    P1P.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_P1P_H
#define radiation_P1P_H

#include "radiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                           Class P1P Declaration
\*---------------------------------------------------------------------------*/

class P1P
:
    public radiationModel
{
    // Private Data

        //- Incident radiation / [W/m^2]
        volScalarField G_;

        //- Total radiative heat flux [W/m^2]
        volScalarField qr_;

        //- Absorption coefficient
        volScalarField a_;

        //- Emission coefficient
        volScalarField e_;

        //- Emission contribution
        volScalarField E_;
        
        //- Particulate equivalent absorption coefficient
        volScalarField ap_;

        //- Equivalent emission of the particles (for calculation here Ep is Ep*pi)
        volScalarField Ep_;
        
        //- G equation left terms first one
        volScalarField GLST1_;


public:

    //- Runtime type information
    TypeName("P1P");


    // Constructors

        //- Construct from components
        P1P(const volScalarField& T);

        //- Construct from components
        P1P(const dictionary& dict, const volScalarField& T);

        //- Disallow default bitwise copy construction
        P1P(const P1P&) = delete;


    //- Destructor
    virtual ~P1P();


    // Member Functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;
            


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const P1P&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
