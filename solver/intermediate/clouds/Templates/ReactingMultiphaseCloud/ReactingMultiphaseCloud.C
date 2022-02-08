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

#include "ReactingMultiphaseCloud.H"

#include "DevolatilisationModel.H"
#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::setModels()
{
    devolatilisationModel_.reset
    (
        DevolatilisationModel<ReactingMultiphaseCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    surfaceReactionModel_.reset
    (
        SurfaceReactionModel<ReactingMultiphaseCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::cloudReset
(
    ReactingMultiphaseCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    devolatilisationModel_.reset(c.devolatilisationModel_.ptr());
    surfaceReactionModel_.reset(c.surfaceReactionModel_.ptr());

    dMassDevolatilisation_ = c.dMassDevolatilisation_;
    dMassSurfaceReaction_ = c.dMassSurfaceReaction_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingMultiphaseCloud<CloudType>::ReactingMultiphaseCloud
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
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    reactingMultiphaseCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    Uavg_(Uavg),
    rhoavg_(rhoavg),
    Tavg_(Tavg),
    muavg_(muavg),
    kappaavg_(kappaavg),
    cpavg_(cpavg),
    alphaavg_(alphaavg),
    o2avg_(o2avg),
    co2avg_(co2avg),
    h2oavg_(h2oavg),
    h2avg_(h2avg),
    devolatilisationModel_(nullptr),
    surfaceReactionModel_(nullptr),
    dMassDevolatilisation_(0.0),
    dMassSurfaceReaction_(0.0)
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
            this->deleteLostParticles();
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ReactingMultiphaseCloud<CloudType>::ReactingMultiphaseCloud
(
    ReactingMultiphaseCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    reactingMultiphaseCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    Uavg_(c.Uavg_),
    rhoavg_(c.rhoavg_),
    Tavg_(c.Tavg_),
    muavg_(c.muavg_),
    kappaavg_(c.kappaavg_),
    cpavg_(c.cpavg_),
    alphaavg_(c.alphaavg_),
    o2avg_(c.o2avg_),
    co2avg_(c.co2avg_),
    h2oavg_(c.h2oavg_),
    h2avg_(c.h2avg_),
    devolatilisationModel_(c.devolatilisationModel_->clone()),
    surfaceReactionModel_(c.surfaceReactionModel_->clone()),
    dMassDevolatilisation_(c.dMassDevolatilisation_),
    dMassSurfaceReaction_(c.dMassSurfaceReaction_)
{}


template<class CloudType>
Foam::ReactingMultiphaseCloud<CloudType>::ReactingMultiphaseCloud
(
    const fvMesh& mesh,
    const word& name,
    const ReactingMultiphaseCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    reactingMultiphaseCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(),
    Uavg_(c.Uavg_),
    rhoavg_(c.rhoavg_),
    Tavg_(c.Tavg_),
    muavg_(c.muavg_),
    kappaavg_(c.kappaavg_),
    cpavg_(c.cpavg_),
    alphaavg_(c.alphaavg_),
    o2avg_(c.o2avg_),
    co2avg_(c.co2avg_),
    h2oavg_(c.h2oavg_),
    h2avg_(c.h2avg_),
    devolatilisationModel_(nullptr),
    surfaceReactionModel_(nullptr),
    dMassDevolatilisation_(0.0),
    dMassSurfaceReaction_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingMultiphaseCloud<CloudType>::~ReactingMultiphaseCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    label idGas = this->composition().idGas();
    label idLiquid = this->composition().idLiquid();
    label idSolid = this->composition().idSolid();

    parcel.YGas() = this->composition().Y0(idGas);
    parcel.YLiquid() = this->composition().Y0(idLiquid);
    parcel.YSolid() = this->composition().Y0(idSolid);
    
    //initialize layer properties from constant properties. Particle contains all layers!
    parcel.Tb0() = constProps_.Tb00();
    parcel.Tb1() = constProps_.Tb10();
    parcel.Tb2() = constProps_.Tb20();
    parcel.Tb3() = constProps_.Tb30();   
    parcel.rb0() = constProps_.rb00();
    parcel.rb1() = constProps_.rb10();
    parcel.rb2() = constProps_.rb20();
    parcel.rb3() = constProps_.rb30();   
    parcel.Tp0() = constProps_.Tp00();
    parcel.Tp1() = constProps_.Tp10();
    parcel.Tp2() = constProps_.Tp20();
    parcel.Tp3() = constProps_.Tp30();
    if(constProps_.parcelShape() == 1)
    {
        parcel.mp0() = constProps_.rho0()*constant::mathematical::pi*(2.0*pow3(parcel.rb0())-sqr(parcel.rb0())*constProps_.xi0()),
        parcel.mp1() = constProps_.rho0()*constant::mathematical::pi*(2.0*pow3(parcel.rb1()) - 2.0*pow3(parcel.rb0()) - (sqr(parcel.rb1()) - sqr(parcel.rb0()))*constProps_.xi0());
        parcel.mp2() = constProps_.rho0()*constant::mathematical::pi*(2.0*pow3(parcel.rb2()) - 2.0*pow3(parcel.rb1()) - (sqr(parcel.rb2()) - sqr(parcel.rb1()))*constProps_.xi0());
        parcel.mp3() = constProps_.rho0()*constant::mathematical::pi*(2.0*pow3(parcel.rb3()) - 2.0*pow3(parcel.rb2()) - (sqr(parcel.rb3()) - sqr(parcel.rb2()))*constProps_.xi0());
    }
    else
    {
        parcel.mp0() = constProps_.rho0()*(4.0/3.0)*constant::mathematical::pi*pow3(parcel.rb0());
        parcel.mp1() = constProps_.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(parcel.rb1())-pow3(parcel.rb0()));
        parcel.mp2() = constProps_.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(parcel.rb2())-pow3(parcel.rb1()));
        parcel.mp3() = constProps_.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(parcel.rb3())-pow3(parcel.rb2()));
    }
    parcel.flagBoiling() = 1;
    parcel.flagDevo() = 1;
        
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    if (fullyDescribed)
    {
        label idGas = this->composition().idGas();
        label idLiquid = this->composition().idLiquid();
        label idSolid = this->composition().idSolid();

        this->checkSuppliedComposition
        (
            parcel.YGas(),
            this->composition().Y0(idGas),
            "YGas"
        );
        this->checkSuppliedComposition
        (
            parcel.YLiquid(),
            this->composition().Y0(idLiquid),
            "YLiquid"
        );
        this->checkSuppliedComposition
        (
            parcel.YSolid(),
            this->composition().Y0(idSolid),
            "YSolid"
        );
    }
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ReactingMultiphaseCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);
Info<<"ReactingMultiphaseCloud<CloudType>::evolve() called here"<<endl;
        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::autoMap
(
    const mapPolyMesh& mapper
)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::info()
{
    CloudType::info();

    this->devolatilisation().info(Info);
    this->surfaceReaction().info(Info);
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::writeFields() const
{
    if (this->compositionModel_.valid())
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
}


// ************************************************************************* //
