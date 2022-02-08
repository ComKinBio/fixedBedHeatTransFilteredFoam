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
inline Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties(),
    TDevol_(this->dict_, 0.0),
    LDevol_(this->dict_, 0.0),
    hRetentionCoeff_(this->dict_, 0.0),
    coarseGrid_(this->dict_, false),
    alphaMax_(this->dict_, 0.0),
    disableHtc_(this->dict_, false),
    noCombustion_(this->dict_, false),
    heatRatioFlag_(this->dict_, false),
    ppConductionFlag_(this->dict_, false),
    controledCombustion_(this->dict_, false),
    combustionAt_(this->dict_, 0.0),
    surfaceExplicit_(this->dict_, false),
    coarseGridSurfaceCombustion_(this->dict_, false),
    radiationCorrection_(this->dict_, 0.0),
    ppRadiationFlag_(this->dict_, false),
    shrinkageMode_(this->dict_, 1),
    shrinkageFactorAsh_(this->dict_, 0.95),
    Tb00_(this->dict_, 0.0),
    Tb10_(this->dict_, 0.0),
    Tb20_(this->dict_, 0.0),
    Tb30_(this->dict_, 0.0),
    rb00_(this->dict_, 0.0),
    rb10_(this->dict_, 0.0),
    rb20_(this->dict_, 0.0),    
    rb30_(this->dict_, 0.0),
    Tp00_(this->dict_, 0.0),
    Tp10_(this->dict_, 0.0),     
    Tp20_(this->dict_, 0.0),  
    Tp30_(this->dict_, 0.0),
    deltaTime_(this->dict_, 0.0),
    hCoeff_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    TDevol_(cp.TDevol_),
    LDevol_(cp.LDevol_),
    hRetentionCoeff_(cp.hRetentionCoeff_),
    coarseGrid_(cp.coarseGrid_),
    alphaMax_(cp.alphaMax_),
    disableHtc_(cp.disableHtc_),
    noCombustion_(cp.noCombustion_),
    heatRatioFlag_(cp.heatRatioFlag_),
    ppConductionFlag_(cp.ppConductionFlag_),
    controledCombustion_(cp.controledCombustion_),
    combustionAt_(cp.combustionAt_),
    surfaceExplicit_(cp.surfaceExplicit_),
    coarseGridSurfaceCombustion_(cp.coarseGridSurfaceCombustion_),
    radiationCorrection_(cp.radiationCorrection_),
    ppRadiationFlag_(cp.ppRadiationFlag_),
    shrinkageMode_(cp.shrinkageMode_),
    shrinkageFactorAsh_(cp.shrinkageFactorAsh_),
    Tb00_(cp.Tb00_),
    Tb10_(cp.Tb10_),
    Tb20_(cp.Tb20_),
    Tb30_(cp.Tb30_),
    rb00_(cp.rb00_),
    rb10_(cp.rb10_),
    rb20_(cp.rb20_),
    rb30_(cp.rb30_),
    Tp00_(cp.Tp00_),
    Tp10_(cp.Tp10_),
    Tp20_(cp.Tp20_),
    Tp30_(cp.Tp30_),
    deltaTime_(cp.deltaTime_),
    hCoeff_(cp.hCoeff_)  
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    TDevol_(this->dict_, "TDevol"),
    LDevol_(this->dict_, "LDevol"),
    hRetentionCoeff_(this->dict_, "hRetentionCoeff"),
    coarseGrid_(this->dict_, word("coarseGrid")),
    alphaMax_(this->dict_, word("alphaMax")),
    disableHtc_(this->dict_, word("disableHtc")),
    noCombustion_(this->dict_, word("noCombustion")),
    heatRatioFlag_(this->dict_, word("useHeatRatio")),
    ppConductionFlag_(this->dict_, word("usePpConduction")),
    controledCombustion_(this->dict_, word("limitedCombustion")),
    combustionAt_(this->dict_, word("devoDegree")),
    surfaceExplicit_(this->dict_, word("surfaceExplicit")),
    coarseGridSurfaceCombustion_(this->dict_, word("coarseGridSpeciesConcentration")),
    radiationCorrection_(this->dict_, word("radiationCorrection")),
    ppRadiationFlag_(this->dict_, word("usePpradiation")),
    shrinkageMode_(this->dict_, word("shrinkageMode")),
    shrinkageFactorAsh_(this->dict_, word("Srinkgamma")),
    Tb00_(this->dict_, "Tb00"),
    Tb10_(this->dict_, "Tb10"),
    Tb20_(this->dict_, "Tb20"),
    Tb30_(this->dict_, "Tb30"),
    rb00_(this->dict_, "rb00"),
    rb10_(this->dict_, "rb10"),
    rb20_(this->dict_, "rb20"),    
    rb30_(this->dict_, "rb30"),
    Tp00_(this->dict_, "Tp00"),
    Tp10_(this->dict_, "Tp10"),     
    Tp20_(this->dict_, "Tp20"),  
    Tp30_(this->dict_, "Tp30"),
    deltaTime_(this->dict_, "deltaTime"),
    hCoeff_(this->dict_, "hCoeff")
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    Tb0_(0.0),
    Tb1_(0.0),
    Tb2_(0.0),
    Tb3_(0.0),
    rb0_(0.0),
    rb1_(0.0),
    rb2_(0.0),    
    rb3_(0.0),
    Tp0_(0.0),
    Tp1_(0.0),     
    Tp2_(0.0),  
    Tp3_(0.0),
    mp0_(0.0),
    mp1_(0.0),
    mp2_(0.0),
    mp3_(0.0),
    rDry_(0.0),
    rDevo_(0.0),
    rChar_(0.0),
    rComb_(0.0),
    ash_inchar_t_(0.0),
    QDry_(0.0),
    QComb_(0.0),
    flagBoiling_(0),
    flagDevo_(0),
    cumTime_(0.0),
    Qin_(0.0),
    QLayer_(0.0),
    VID1_(0.0),
    VID2_(0.0),
    VID3_(0.0),
    VID4_(0.0),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    Tb0_(0.0),
    Tb1_(0.0),
    Tb2_(0.0),
    Tb3_(0.0),
    rb0_(0.0),
    rb1_(0.0),
    rb2_(0.0),    
    rb3_(0.0),
    Tp0_(0.0),
    Tp1_(0.0),     
    Tp2_(0.0),  
    Tp3_(0.0),
    mp0_(0.0),
    mp1_(0.0),
    mp2_(0.0),
    mp3_(0.0),
    rDry_(0.0),
    rDevo_(0.0),
    rChar_(0.0),
    rComb_(0.0),
    ash_inchar_t_(0.0),
    QDry_(0.0),
    QComb_(0.0),
    flagBoiling_(0),
    flagDevo_(0),
    cumTime_(0.0),
    Qin_(0.0),
    QLayer_(0.0),
    VID1_(0.0),
    VID2_(0.0),
    VID3_(0.0),
    VID4_(0.0),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
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
    const scalarField& Y0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
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
        Y0,
        constProps
    ),
    Tb0_(constProps.Tb00()),
    Tb1_(constProps.Tb10()),
    Tb2_(constProps.Tb20()),
    Tb3_(constProps.Tb30()),
    rb0_(constProps.rb00()),
    rb1_(constProps.rb10()),
    rb2_(constProps.rb20()),    
    rb3_(constProps.rb30()),
    Tp0_(constProps.Tp00()),
    Tp1_(constProps.Tp10()),     
    Tp2_(constProps.Tp20()),  
    Tp3_(constProps.Tp30()),
    mp0_(constProps.rho0()*(4.0/3.0)*constant::mathematical::pi*pow3(rb0_)),
    mp1_(constProps.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(rb1_)-pow3(rb0_))),
    mp2_(constProps.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(rb2_)-pow3(rb1_))),
    mp3_(constProps.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(rb3_)-pow3(rb2_))),
    rDry_(0.0),
    rDevo_(0.0),
    rChar_(0.0),
    rComb_(0.0),
    ash_inchar_t_(0.0),
    QDry_(0.0),
    QComb_(0.0),
    flagBoiling_(1),
    flagDevo_(1),
    cumTime_(0.0),
    Qin_(0.0),
    QLayer_(0.0), 
    VID1_(0.0),
    VID2_(0.0),
    VID3_(0.0),
    VID4_(0.0),
    YGas_(YGas0),
    YLiquid_(YLiquid0),
    YSolid_(YSolid0),
    canCombust_(0)
{
    if(constProps.parcelShape() == 1)
    {
        mp0_ = constProps.rho0()*constant::mathematical::pi*(2*pow3(rb0_)-sqr(rb0_)*constProps.xi0()),
        mp1_ = constProps.rho0()*constant::mathematical::pi*(2.0*pow3(rb1_) - 2.0*pow3(rb0_) - (sqr(rb1_) - sqr(rb0_))*constProps.xi0());
        mp2_ = constProps.rho0()*constant::mathematical::pi*(2.0*pow3(rb2_) - 2.0*pow3(rb1_) - (sqr(rb2_) - sqr(rb1_))*constProps.xi0());
        mp3_ = constProps.rho0()*constant::mathematical::pi*(2.0*pow3(rb3_) - 2.0*pow3(rb2_) - (sqr(rb3_) - sqr(rb2_))*constProps.xi0());
    }
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::TDevol() const
{
    return TDevol_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::LDevol() const
{
    return LDevol_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
hRetentionCoeff() const
{
    scalar value = hRetentionCoeff_.value();

    if ((value < 0) || (value > 1))
    {
        FatalErrorInFunction
            << "hRetentionCoeff must be in the range 0 to 1" << nl
            << exit(FatalError) << endl;
    }

    return value;
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Tb00() const
{
    return Tb00_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Tb10() const
{
    return Tb10_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Tb20() const
{
    return Tb20_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Tb30() const
{
    return Tb30_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::rb00() const
{
    return rb00_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::rb10() const
{
    return rb10_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::rb20() const
{
    return rb20_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::rb30() const
{
    return rb30_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Tp00() const
{
    return Tp00_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Tp10() const
{
    return Tp10_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Tp20() const
{
    return Tp20_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Tp30() const
{
    return Tp30_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::deltaTime() const
{
    return deltaTime_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::hCoeff() const
{
    return hCoeff_.value();
}

template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::coarseGrid() const
{
    return coarseGrid_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::alphaMax() const
{
    return alphaMax_.value();
}

template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::disableHtc() const
{
    return disableHtc_.value();
}

template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::noCombustion() const
{
    return noCombustion_.value();
}

template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::heatRatioFlag() const
{
    return heatRatioFlag_.value();
}


template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::ppConductionFlag() const
{
    return ppConductionFlag_.value();
}

template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::controledCombustion() const
{
    return controledCombustion_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::combustionAt() const
{
    return combustionAt_.value();
}

template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::surfaceExplicit() const
{
    return surfaceExplicit_.value();
}

template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::coarseGridSurfaceCombustion() const
{
    return coarseGridSurfaceCombustion_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::radiationCorrection() const
{
    return radiationCorrection_.value();
}

template<class ParcelType>
inline bool
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::ppRadiationFlag() const
{
    return ppRadiationFlag_.value();
}

template<class ParcelType>
inline Foam::label
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::shrinkageMode() const
{
    return shrinkageMode_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::shrinkageFactorAsh() const
{
    return shrinkageFactorAsh_.value();
}

// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YGas() const
{
    return YGas_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YLiquid() const
{
    return YLiquid_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YSolid() const
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::label
Foam::ReactingMultiphaseParcel<ParcelType>::canCombust() const
{
    return canCombust_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Tb0() const
{
    return Tb0_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Tb1() const
{
    return Tb1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Tb2() const
{
    return Tb2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Tb3() const
{
    return Tb3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::rb0() const
{
    return rb0_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::rb1() const
{
    return rb1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::rb2() const
{
    return rb2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::rb3() const
{
    return rb3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Tp0() const
{
    return Tp0_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Tp1() const
{
    return Tp1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Tp2() const
{
    return Tp2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Tp3() const
{
    return Tp3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::mp0() const
{
    return mp0_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::mp1() const
{
    return mp1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::mp2() const
{
    return mp2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::mp3() const
{
    return mp3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::rDry() const
{
    return rDry_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::rDevo() const
{
    return rDevo_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::rChar() const
{
    return rChar_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::rComb() const
{
    return rComb_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::ash_inchar_t() const
{
    return ash_inchar_t_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::QDry() const
{
    return QDry_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::QComb() const
{
    return QComb_;
}

template<class ParcelType>
inline Foam::label
Foam::ReactingMultiphaseParcel<ParcelType>::flagBoiling() const
{
    return flagBoiling_;
}

template<class ParcelType>
inline Foam::label
Foam::ReactingMultiphaseParcel<ParcelType>::flagDevo() const
{
    return flagDevo_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::cumTime() const
{
    return cumTime_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::Qin() const
{
    return Qin_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::QLayer() const
{
    return QLayer_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::VID1() const
{
    return VID1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::VID2() const
{
    return VID2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::VID3() const
{
    return VID3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::VID4() const
{
    return VID4_;
}

template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YGas()
{
    return YGas_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YLiquid()
{
    return YLiquid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YSolid()
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::label& Foam::ReactingMultiphaseParcel<ParcelType>::canCombust()
{
    return canCombust_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Tb0()
{
    return Tb0_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Tb1()
{
    return Tb1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Tb2()
{
    return Tb2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Tb3()
{
    return Tb3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::rb0()
{
    return rb0_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::rb1()
{
    return rb1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::rb2()
{
    return rb2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::rb3()
{
    return rb3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Tp0()
{
    return Tp0_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Tp1()
{
    return Tp1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Tp2()
{
    return Tp2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Tp3()
{
    return Tp3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::mp0()
{
    return mp0_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::mp1()
{
    return mp1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::mp2()
{
    return mp2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::mp3()
{
    return mp3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::rDry()
{
    return rDry_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::rDevo()
{
    return rDevo_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::rChar()
{
    return rChar_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::rComb()
{
    return rComb_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::ash_inchar_t()
{
    return ash_inchar_t_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::QDry()
{
    return QDry_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::QComb()
{
    return QComb_;
}

template<class ParcelType>
inline Foam::label&
Foam::ReactingMultiphaseParcel<ParcelType>::flagBoiling()
{
    return flagBoiling_;
}

template<class ParcelType>
inline Foam::label&
Foam::ReactingMultiphaseParcel<ParcelType>::flagDevo()
{
    return flagDevo_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::cumTime()
{
    return cumTime_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::Qin()
{
    return Qin_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::QLayer()
{
    return QLayer_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::VID1()
{
    return VID1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::VID2()
{
    return VID2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::VID3()
{
    return VID3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseParcel<ParcelType>::VID4()
{
    return VID4_;
}

// ************************************************************************* //
