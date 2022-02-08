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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    T0_(this->dict_, 0.0),
    TMin_(this->dict_, 0.0),
    TMax_(this->dict_, vGreat),
    Cp0_(this->dict_, 0.0),
    kp0_(this->dict_, 0.0),
    epsilon0_(this->dict_, 0.0),
    f0_(this->dict_, 0.0),
    ppConductionFlag_(this->dict_, false),
    ppRadiationFlag_(this->dict_, false)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    T0_(cp.T0_),
    TMin_(cp.TMin_),
    TMax_(cp.TMax_),
    Cp0_(cp.Cp0_),
    kp0_(cp.kp0_),
    epsilon0_(cp.epsilon0_),
    f0_(cp.f0_),
    ppConductionFlag_(cp.ppConductionFlag_),
    ppRadiationFlag_(cp.ppConductionFlag_)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    T0_(this->dict_, "T0"),
    TMin_(this->dict_, "TMin", 200.0),
    TMax_(this->dict_, "TMax", 5000.0),
    Cp0_(this->dict_, "Cp0"),
    kp0_(this->dict_, "kp0"),
    epsilon0_(this->dict_, "epsilon0"),
    f0_(this->dict_, "f0"),
    ppConductionFlag_(this->dict_, "usePpConduction", false),
    ppRadiationFlag_(this->dict_, "usePpradiation", false)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    T_(0.0),
    Cp_(0.0),
    kp_(0.0),
    dT_(0.0),
    ppT4sum_(0.0),
    neighbourNum_(0.0),
    neighbourMax_(0.0),
    neighbourMin_(0.0),
    dBeforeDEM_(0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    T_(0.0),
    Cp_(0.0),
    kp_(0.0),
    dT_(0.0),
    ppT4sum_(0.0),
    neighbourNum_(0.0),
    neighbourMax_(0.0),
    neighbourMin_(0.0),
    dBeforeDEM_(0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    T_(constProps.T0()),
    Cp_(constProps.Cp0()),
    kp_(constProps.kp0()),
    dT_(0.0),
    ppT4sum_(0.0),
    neighbourNum_(0.0),
    neighbourMax_(0.0),
    neighbourMin_(0.0),
    dBeforeDEM_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::T0() const
{
    return T0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::TMin() const
{
    return TMin_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::TMax() const
{
    return TMax_.value();
}


template<class ParcelType>
inline void
Foam::ThermoParcel<ParcelType>::constantProperties::setTMax(const scalar TMax)
{
    TMax_.setValue(TMax);
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::Cp0() const
{
    return Cp0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::kp0() const
{
    return kp0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::epsilon0() const
{
    return epsilon0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::f0() const
{
    return f0_.value();
}


template<class ParcelType>
inline bool
Foam::ThermoParcel<ParcelType>::constantProperties::ppConductionFlag() const
{
    return ppConductionFlag_.value();
}


template<class ParcelType>
inline bool
Foam::ThermoParcel<ParcelType>::constantProperties::ppRadiationFlag() const
{
    return ppRadiationFlag_.value();
}

// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::kp() const
{
    return kp_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::hs() const
{
    return Cp_*(T_ - 298.15);
}

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::dT() const
{
    return dT_;
}

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::ppT4sum() const
{
    return ppT4sum_;
}

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::neighbourNum() const
{
    return neighbourNum_;
}

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::neighbourMax() const
{
    return neighbourMax_;
}

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::neighbourMin() const
{
    return  neighbourMin_;
}

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::dBeforeDEM() const
{
    return  dBeforeDEM_;
}

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::kpp(const scalar T) const
{
    return kp();
}

template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::T()
{
    return T_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::dT()
{
    return dT_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::ppT4sum()
{
    return ppT4sum_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::neighbourNum()
{
    return neighbourNum_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::neighbourMax()
{
    return neighbourMax_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::neighbourMin()
{
    return  neighbourMin_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::dBeforeDEM()
{
    return  dBeforeDEM_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::Cp()
{
    return Cp_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::kp()
{
    return kp_;
}


// ************************************************************************* //
