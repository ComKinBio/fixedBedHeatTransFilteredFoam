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
    Foam::ReactingMultiphaseCloud

Description
    Templated base class for multiphase reacting cloud

    - Adds to reacting cloud
      - multiphase composition
      - devolatilisatsion
      - surface reactions

SourceFiles
    ReactingMultiphaseCloudI.H
    ReactingMultiphaseCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseCloud_H
#define ReactingMultiphaseCloud_H

#include "reactingMultiphaseCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class DevolatilisationModel;

template<class CloudType>
class SurfaceReactionModel;

/*---------------------------------------------------------------------------*\
                   Class ReactingMultiphaseCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReactingMultiphaseCloud
:
    public CloudType,
    public reactingMultiphaseCloud
{
public:

    // Public Typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef ReactingMultiphaseCloud<CloudType> reactingMultiphaseCloudType;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<ReactingMultiphaseCloud<CloudType>> cloudCopyPtr_;


protected:

    // Protected data

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;

        //- goarse grid average velocity [m/s]
        const volVectorField& Uavg_;
        
        //- goarse grid average density 
        const volScalarField& rhoavg_;
        
        //- goarse grid average Temperature [K]
        const volScalarField& Tavg_;
        
        //- goarse grid average viscosity [m/s]
        const volScalarField& muavg_;
        
        //- goarse grid average thermal conductivity [m/s]
        const volScalarField& kappaavg_;
        
        //- goarse grid average heat capacity [m/s]
        const volScalarField& cpavg_;
        
        //- goarse grid average solid phase fraction 
        const volScalarField& alphaavg_;
        
        //- goarse grid average gas phase concentration mass fraction O2
        const volScalarField& o2avg_;
        
        //- goarse grid average gas phase concentration mass fraction CO2
        const volScalarField& co2avg_;
        
        //- goarse grid average gas phase concentration mass fraction H2O
        const volScalarField& h2oavg_;
        
        //- goarse grid average gas phase concentration mass fraction H2
        const volScalarField& h2avg_;

        // References to the cloud sub-models

            //- Devolatilisation model
            autoPtr
            <
                DevolatilisationModel<ReactingMultiphaseCloud<CloudType>>
            >
            devolatilisationModel_;

            //- Surface reaction model
            autoPtr
            <
                SurfaceReactionModel<ReactingMultiphaseCloud<CloudType>>
            >
            surfaceReactionModel_;


        // Check

            //- Total mass transferred to continuous phase via devolatilisation
            scalar dMassDevolatilisation_;

            //- Total mass transferred to continuous phase via surface
            //  reactions
            scalar dMassSurfaceReaction_;


    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(ReactingMultiphaseCloud<CloudType>& c);


public:

    // Constructors

//         //- Construct given carrier gas fields
//         ReactingMultiphaseCloud
//         (
//             const word& cloudName,
//             const volScalarField& rho,
//             const volVectorField& U,
//             const dimensionedVector& g,
//             const SLGThermo& thermo,
//             bool readFields = true
//         );
        
        //- Construct given carrier gas fields
        ReactingMultiphaseCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            const volVectorField& Uavg,
            const volScalarField& rhoavg,
            const volScalarField& Tavg,
            const volScalarField& muavg,
            const volScalarField& kappaavg,
            const volScalarField& cpavg,
            const volScalarField& alphaavg,
            const volScalarField& o2avg,
            const volScalarField& co2avg,
            const volScalarField& h2oavg,
            const volScalarField& h2avg,
            bool readFields = true
        );


        //- Copy constructor with new name
        ReactingMultiphaseCloud
        (
            ReactingMultiphaseCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        ReactingMultiphaseCloud
        (
            const fvMesh& mesh,
            const word& name,
            const ReactingMultiphaseCloud<CloudType>& c
        );

        //- Disallow default bitwise copy construction
        ReactingMultiphaseCloud(const ReactingMultiphaseCloud&) = delete;

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ReactingMultiphaseCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ReactingMultiphaseCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~ReactingMultiphaseCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const ReactingMultiphaseCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();
            
            //- Return carrier gas velocity
            inline const volVectorField& Uavg() const;

            //- Return carrier gas density
            inline const volScalarField& rhoavg() const;

            //- Return carrier gas dynamic viscosity
            inline const volScalarField& muavg() const;
            
            //- Return const access to the carrier average temperature field
            inline const volScalarField& Tavg() const;
            
            //- Return const access to the carrier average thermal conductivity field
            inline const volScalarField& kappaavg() const;
            
            //- Return const access to the carrier average heat capacity field
            inline const volScalarField& cpavg() const;
            
            //- Return const access to the carrier average temperature field
            inline const volScalarField& alphaavg() const;
            
            //- Return const access to the carrier average temperature field
            inline const volScalarField& o2avg() const;
            
            //- Return const access to the carrier average thermal conductivity field
            inline const volScalarField& h2oavg() const;
            
            //- Return const access to the carrier average heat capacity field
            inline const volScalarField& co2avg() const;
            
            //- Return const access to the carrier average temperature field
            inline const volScalarField& h2avg() const;

            // Sub-models

                //- Return const access to devolatilisation model
                inline const DevolatilisationModel
                <
                    ReactingMultiphaseCloud<CloudType>
                >&
                devolatilisation() const;

                //- Return reference to devolatilisation model
                inline DevolatilisationModel
                <
                    ReactingMultiphaseCloud<CloudType>
                >&
                devolatilisation();

                //- Return const access to reacting surface reaction model
                inline const SurfaceReactionModel
                <
                    ReactingMultiphaseCloud<CloudType>
                >&
                surfaceReaction() const;

                //- Return reference to reacting surface reaction model
                inline SurfaceReactionModel
                <
                    ReactingMultiphaseCloud<CloudType>
                >&
                surfaceReaction();


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Evolve the cloud
            void evolve();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();

            //- Write the field data for the cloud
            virtual void writeFields() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ReactingMultiphaseCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingMultiphaseCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingMultiphaseCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
