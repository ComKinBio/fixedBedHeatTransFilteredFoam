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

#include "PairSpringSliderDashpotHeat.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::PairSpringSliderDashpotHeat<CloudType>::findMinMaxProperties
(
    scalar& RMin,
    scalar& rhoMax,
    scalar& UMagMax
) const
{
    RMin = vGreat;
    rhoMax = -vGreat;
    UMagMax = -vGreat;

    forAllConstIter(typename CloudType, this->owner(), iter)
    {
        const typename CloudType::parcelType& p = iter();

        // Finding minimum diameter to avoid excessive arithmetic

        scalar dEff = p.d();

        if (useEquivalentSize_)
        {
            dEff *= cbrt(p.nParticle()*volumeFactor_);
        }

        RMin = min(dEff, RMin);

        rhoMax = max(p.rho(), rhoMax);

        UMagMax = max
        (
            mag(p.U()) + mag(p.omega())*dEff/2,
            UMagMax
        );
    }

    // Transform the minimum diameter into minimum radius
    //     rMin = dMin/2
    // then rMin into minimum R,
    //     1/RMin = 1/rMin + 1/rMin,
    //     RMin = rMin/2 = dMin/4

    RMin /= 4.0;

    // Multiply by two to create the worst-case relative velocity

    UMagMax = 2*UMagMax;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairSpringSliderDashpotHeat<CloudType>::PairSpringSliderDashpotHeat
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PairModel<CloudType>(dict, cloud, typeName),
    Estar_(),
    Gstar_(),
    alpha_(readScalar(this->coeffDict().lookup("alpha"))),
    b_(readScalar(this->coeffDict().lookup("b"))),
    mu_(readScalar(this->coeffDict().lookup("mu"))),
    cohesionEnergyDensity_
    (
        readScalar(this->coeffDict().lookup("cohesionEnergyDensity"))
    ),
    cohesion_(false),
    collisionResolutionSteps_
    (
        readScalar(this->coeffDict().lookup("collisionResolutionSteps"))
    ),
    volumeFactor_(1.0),
    useEquivalentSize_(Switch(this->coeffDict().lookup("useEquivalentSize")))
{
    if (useEquivalentSize_)
    {
        volumeFactor_ = readScalar(this->coeffDict().lookup("volumeFactor"));
    }

    scalar nu = this->owner().constProps().poissonsRatio();

    scalar E = this->owner().constProps().youngsModulus();

    Estar_ = E/(2.0*(1.0 - sqr(nu)));

    scalar G = E/(2.0*(1.0 + nu));

    Gstar_ = G/(2.0*(2.0 - nu));

    cohesion_ = (mag(cohesionEnergyDensity_) > vSmall);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairSpringSliderDashpotHeat<CloudType>::~PairSpringSliderDashpotHeat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PairSpringSliderDashpotHeat<CloudType>::controlsTimestep() const
{
    return true;
}


template<class CloudType>
Foam::label Foam::PairSpringSliderDashpotHeat<CloudType>::nSubCycles() const
{
    if (!(this->owner().size()))
    {
        return 1;
    }

    scalar RMin;
    scalar rhoMax;
    scalar UMagMax;

    findMinMaxProperties(RMin, rhoMax, UMagMax);
    
    //TODO need to be checked
    if (UMagMax<=vSmall)
    {
        UMagMax = 0.01;
    }
    
    // Note:  pi^(7/5)*(5/4)^(2/5) = 5.429675
    scalar minCollisionDeltaT =
        5.429675
       *RMin
       *pow(rhoMax/(Estar_*sqrt(UMagMax) + vSmall), 0.4)
       /collisionResolutionSteps_;

    return ceil(this->owner().time().deltaTValue()/minCollisionDeltaT);
}


template<class CloudType>
void Foam::PairSpringSliderDashpotHeat<CloudType>::evaluatePair
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    vector r_AB = (pA.position() - pB.position());

    scalar dAEff = pA.d();

    if (useEquivalentSize_)
    {
        dAEff *= cbrt(pA.nParticle()*volumeFactor_);
    }

    scalar dBEff = pB.d();

    if (useEquivalentSize_)
    {
        dBEff *= cbrt(pB.nParticle()*volumeFactor_);
    }

    scalar r_AB_mag = mag(r_AB);

    scalar normalOverlapMag = 0.5*(dAEff + dBEff) - r_AB_mag;

    if (normalOverlapMag > 0)
    {
        // Particles in collision

        vector rHat_AB = r_AB/(r_AB_mag + vSmall);

        vector U_AB = pA.U() - pB.U();

        // Effective radius
        scalar R = 0.5*dAEff*dBEff/(dAEff + dBEff);

        // Effective mass
        scalar M = pA.mass()*pB.mass()/(pA.mass() + pB.mass());

        scalar kN = (4.0/3.0)*sqrt(R)*Estar_;

        scalar etaN = alpha_*sqrt(M*kN)*pow025(normalOverlapMag);

        // Normal force
        vector fN_AB =
            rHat_AB
           *(kN*pow(normalOverlapMag, b_) - etaN*(U_AB & rHat_AB));

        // Cohesion force, energy density multiplied by the area of
        // particle-particle overlap
        if (cohesion_)
        {
            fN_AB +=
                -cohesionEnergyDensity_
                *overlapArea(dAEff/2.0, dBEff/2.0, r_AB_mag)
                *rHat_AB;
        }

        pA.f() += fN_AB;
        pB.f() += -fN_AB;

        vector USlip_AB =
            U_AB - (U_AB & rHat_AB)*rHat_AB
          - ((dAEff/2*pA.omega() + dBEff/2*pB.omega()) ^ rHat_AB);

        scalar deltaT = this->owner().mesh().time().deltaTValue();

        vector& tangentialOverlap_AB =
            pA.collisionRecords().matchPairRecord
            (
                pB.origProc(),
                pB.origId()
            ).collisionData();

        vector& tangentialOverlap_BA =
            pB.collisionRecords().matchPairRecord
            (
                pA.origProc(),
                pA.origId()
            ).collisionData();

        vector deltaTangentialOverlap_AB = USlip_AB*deltaT;

        tangentialOverlap_AB += deltaTangentialOverlap_AB;
        tangentialOverlap_BA += -deltaTangentialOverlap_AB;

        scalar tangentialOverlapMag = mag(tangentialOverlap_AB);

        if (tangentialOverlapMag > vSmall)
        {
            scalar kT = 8.0*sqrt(R*normalOverlapMag)*Gstar_;

            scalar etaT = etaN;

            // Tangential force
            vector fT_AB;

            if (kT*tangentialOverlapMag > mu_*mag(fN_AB))
            {
                // Tangential force greater than sliding friction,
                // particle slips

                fT_AB = -mu_*mag(fN_AB)*USlip_AB/mag(USlip_AB);

                tangentialOverlap_AB = Zero;
                tangentialOverlap_BA = Zero;
            }
            else
            {
                fT_AB = - kT*tangentialOverlap_AB - etaT*USlip_AB;
            }

            pA.f() += fT_AB;
            pB.f() += -fT_AB;

            pA.torque() += (dAEff/2*-rHat_AB) ^ fT_AB;
            pB.torque() += (dBEff/2*rHat_AB) ^ -fT_AB;
        }
    }
}

template<class CloudType>
void Foam::PairSpringSliderDashpotHeat<CloudType>::evaluatePairHeat
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    vector r_AB = (pA.position() - pB.position());
    
    scalar r_AB_mag = mag(r_AB);
    
    scalar dAEff = pA.dBeforeDEM();

    if (useEquivalentSize_)
    {
        dAEff *= cbrt(pA.nParticle()*volumeFactor_);
    }

    scalar dBEff = pB.dBeforeDEM();

    if (useEquivalentSize_)
    {
        dBEff *= cbrt(pB.nParticle()*volumeFactor_);
    }   
   
    scalar normalOverlapMag = 0.5*(dAEff + dBEff) - r_AB_mag;

    if (normalOverlapMag > 0)
    {
        // Particles in collision
     
        //update partice contacted neighbour number
        pA.neighbourNum() += 1.0;
        pB.neighbourNum() += 1.0;
       
        //update partice contacted T4
        pA.ppT4sum() += Foam::pow4(pB.T());
        pB.ppT4sum() += Foam::pow4(pA.T());
       
        //update particle neighbour Max and Min T
        //this code dont have real functions just for analysis
        if(pB.T() > pA.neighbourMax())
        {
            pA.neighbourMax() = pB.T();
        }
        if(pB.T() < pA.neighbourMin())
        {
            pA.neighbourMin() = pB.T();
        }
        
        if(pA.T() > pB.neighbourMax())
        {
            pB.neighbourMax() = pA.T();
        }
        if(pA.T() < pB.neighbourMin())
        {
            pB.neighbourMin() = pA.T();
        }
      
        //calculate heat conduction
        scalar Q_pBtopA = 0.0;
        
        scalar Q_PtoFtoP = 0; //HOLD
            
            //temperature difference between particles
            scalar dT_AB = pB.T() - pA.T();          
            //thermal conductivity (W m-1 K-1)
            scalar k_p_A = pA.kpp(pA.T());
            scalar k_p_B = pB.kpp(pB.T());          
            
            //Heron's formula to calculate radius of the contact area of two contacting particles
            const scalar HF_s = (r_AB_mag + 0.5*(dAEff + dBEff))/2.;
            const scalar HF_A = sqrt(mag(HF_s*(HF_s - r_AB_mag)*(HF_s - dBEff/2.)*(HF_s - dAEff/2.)));
            const scalar r_c = 2.0*HF_A/r_AB_mag;

            //Batchelor and O'Brian 1977.  Zhou and Yu AIChE 2009, eq 16
            Q_pBtopA = (4 * r_c * dT_AB) / (1/k_p_A + 1/k_p_B);
            

            scalar Q_conduct = Q_pBtopA + Q_PtoFtoP;

            //update dT_
            // pp.dT_ = pp.dT_ + Q_conduct / (pp.mass() * pp.Cp());
            pA.dT() += Q_conduct;
            pB.dT() += -Q_conduct;
            
    }
}


// ************************************************************************* //
// if(std::isnan(pA.dT()) || std::isnan(pB.dT()))
// {
//     std::cout<<"processor: "<<Pstream::myProcNo()<<' '<<"particle dT is nan"<<'\n';
//     std::cout<<"partice A position is : "<<"("<<pA.position().x()<<","<<pA.position().y() <<","<<pA.position().z()<<")"<<' '<<"diameter is : "<<dAEff <<' '<<"temperature is : "<<pA.T() <<' '<<"Tp3 is : "<<pA.Tp3() <<' '<<"k_p_A is : "<<k_p_A <<' '<<'\n';
//     std::cout<<"partice B position is : "<<"("<<pB.position().x()<<","<<pB.position().y() <<","<<pB.position().z()<<")"<<' '<<"diameter is : "<<dBEff <<' '<<"temperature is : "<<pB.T() <<' '<<"Tp3 is : "<<pB.Tp3() <<' '<<"k_p_B is : "<<k_p_B <<' '<<'\n';
//     std::cout<<"r_AB_mag is : "<<r_AB_mag <<'\n';
// }
 
