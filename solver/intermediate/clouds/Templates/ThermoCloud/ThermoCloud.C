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

#include "ThermoCloud.H"
#include "ThermoParcel.H"

#include "HeatTransferModel.H"

#include "CollisionModel.H"
#include "NoCollision.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setModels()
{
    heatTransferModel_.reset
    (
        HeatTransferModel<ThermoCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    TIntegrator_.reset
    (
        integrationScheme::New
        (
            "T",
            this->solution().integrationSchemes()
        ).ptr()
    );

    if (this->solution().coupled())
    {
        this->subModelProperties().lookup("radiation") >> radiation_;
    }

    if (radiation_)
    {
        radAreaP_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + ":radAreaP",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimArea, 0)
            )
        );

        radT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + ":radT4",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(pow4(dimTemperature), 0)
            )
        );

        radAreaPT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + ":radAreaPT4",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(sqr(dimLength)*pow4(dimTemperature), 0)
            )
        );
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::cloudReset(ThermoCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    heatTransferModel_.reset(c.heatTransferModel_.ptr());
    TIntegrator_.reset(c.TIntegrator_.ptr());

    radiation_ = c.radiation_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        rho,
        U,
        thermo.thermo().mu(),
        g,
        false
    ),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    thermo_(thermo),
    T_(thermo.thermo().T()),
    p_(thermo.thermo().p()),
    Tavg_(thermo.thermo().T()),
    kappaavg_(thermo.thermo().kappa()),
    cpavg_(thermo.thermo().Cp()),
    heatTransferModel_(nullptr),
    TIntegrator_(nullptr),
    radiation_(false),
    radAreaP_(nullptr),
    radT4_(nullptr),
    radAreaPT4_(nullptr),
    hsTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy, 0)
        )
    ),
    hsCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTemperature, 0)
        )
    ),
    aDispAvg_
    (
        IOobject
        (
            this->name() + "aDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    ),
    eDispAvg_
    (
        IOobject
        (
            this->name() + "eDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    ),
    EDispAvg_
    (
        IOobject
        (
            this->name() + "EDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    ),
    sigmapAvg_
    (
        IOobject
        (
            this->name() +"sigmapAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    )
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}

template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    const volVectorField& Uavg,
    const volScalarField& rhoavg,
    const volScalarField& muavg,
    const volScalarField& alphaavg,
    const volScalarField& Tavg,
    const volScalarField& kappaavg,
    const volScalarField& cpavg,
    sourceFilter* sourceFilterModel,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        rho,
        U,
        thermo.thermo().mu(),
        g,
        Uavg,
        rhoavg,
        muavg,
        alphaavg,
        sourceFilterModel,
        false
    ),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    thermo_(thermo),
    T_(thermo.thermo().T()),
    p_(thermo.thermo().p()),
    Tavg_(Tavg),
    kappaavg_(kappaavg),
    cpavg_(cpavg),
    heatTransferModel_(nullptr),
    TIntegrator_(nullptr),
    radiation_(false),
    radAreaP_(nullptr),
    radT4_(nullptr),
    radAreaPT4_(nullptr),
    hsTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy, 0)
        )
    ),
    hsCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTemperature, 0)
        )
    ),
    aDispAvg_
    (
        IOobject
        (
            this->name() + "aDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    ),
    eDispAvg_
    (
        IOobject
        (
            this->name() + "eDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    ),
    EDispAvg_
    (
        IOobject
        (
            this->name() + "EDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    ),
    sigmapAvg_
    (
        IOobject
        (
            this->name() +"sigmapAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    )
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    ThermoCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    Tavg_(c.Tavg()),
    kappaavg_(c.kappaavg()),
    cpavg_(c.cpavg()),
    heatTransferModel_(c.heatTransferModel_->clone()),
    TIntegrator_(c.TIntegrator_->clone()),
    radiation_(c.radiation_),
    radAreaP_(nullptr),
    radT4_(nullptr),
    radAreaPT4_(nullptr),
    hsTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsTrans()
        )
    ),
    hsCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsCoeff()
        )
    ),
    aDispAvg_
    (
        IOobject
        (
            this->name() + "aDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        c.aDispAvg()
    ),
    eDispAvg_
    (
        IOobject
        (
            this->name() + "eDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        c.eDispAvg()
    ),
    EDispAvg_
    (
        IOobject
        (
            this->name() + "EDispAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        c.EDispAvg()
    ),
    sigmapAvg_
    (
        IOobject
        (
            this->name() +"sigmapAvg",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        c.sigmapAvg()
    )
{
    if (radiation_)
    {
        radAreaP_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + ":radAreaP",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.radAreaP()
            )
        );

        radT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + ":radT4",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.radT4()
            )
        );

        radAreaPT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + ":radAreaPT4",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.radAreaPT4()
            )
        );
    }
}


template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const fvMesh& mesh,
    const word& name,
    const ThermoCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(),
    thermo_(c.thermo()),
    T_(c.T()),
    p_(c.p()),
    Tavg_(c.Tavg()),
    kappaavg_(c.kappaavg()),
    cpavg_(c.cpavg()),
    heatTransferModel_(nullptr),
    TIntegrator_(nullptr),
    radiation_(false),
    radAreaP_(nullptr),
    radT4_(nullptr),
    radAreaPT4_(nullptr),
    hsTrans_(nullptr),
    hsCoeff_(nullptr),
    aDispAvg_(c.aDispAvg()),
    eDispAvg_(c.eDispAvg()),
    EDispAvg_(c.EDispAvg()),
    sigmapAvg_(c.sigmapAvg())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoCloud<CloudType>::~ThermoCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.T() = constProps_.T0();
    parcel.Cp() = constProps_.Cp0();
    parcel.kp() = constProps_.kp0();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ThermoCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::filerSourceTerms()
{
    CloudType::filerSourceTerms();
    if (this->sourceFilterModel().getHeatDiffusionFlag())
    {
        hsTrans() = this->sourceFilterModel().diffusion(hsTrans(),"heat");
        hsCoeff() = this->sourceFilterModel().diffusion(hsCoeff(),"heat");
    }
   
    if (radiation_)
    {
        if (this->sourceFilterModel().getRadiaDiffusionFlag())
        {
            radAreaP() = this->sourceFilterModel().diffusion(radAreaP(),"radia");
            radT4() = this->sourceFilterModel().diffusion(radT4(),"radia");
            radAreaPT4() = this->sourceFilterModel().diffusion(radAreaPT4(),"radia");
        }
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    hsTrans_->field() = 0.0;
    hsCoeff_->field() = 0.0;

    if (radiation_)
    {
        radAreaP_->field() = 0.0;
        radT4_->field() = 0.0;
        radAreaPT4_->field() = 0.0;
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::relaxSources
(
    const ThermoCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    this->relax(hsTrans_(), cloudOldTime.hsTrans(), "h");
    this->relax(hsCoeff_(), cloudOldTime.hsCoeff(), "h");

    if (radiation_)
    {
        this->relax(radAreaP_(), cloudOldTime.radAreaP(), "radiation");
        this->relax(radT4_(), cloudOldTime.radT4(), "radiation");
        this->relax(radAreaPT4_(), cloudOldTime.radAreaPT4(), "radiation");
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(hsTrans_(), "h");
    this->scale(hsCoeff_(), "h");

    if (radiation_)
    {
        this->scale(radAreaP_(), "radiation");
        this->scale(radT4_(), "radiation");
        this->scale(radAreaPT4_(), "radiation");
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::preEvolve()
{
    CloudType::preEvolve();

    this->pAmbient() = thermo_.thermo().p().average().value();//error happens here!!!!!!!!!!!!!TODO
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::info()
{
    CloudType::info();

    Info<< "    Temperature min/max             = " << Tmin() << ", " << Tmax()
        << endl;
}


// ************************************************************************* //
