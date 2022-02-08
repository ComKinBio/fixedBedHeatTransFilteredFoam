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

#include "ReactingMultiphaseParcel.H"
#include "mathematicalConstants.H"


#define CP_GAS 1.18698e3

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::GAS(0);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::LIQ(1);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::SLD(2);


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::CpEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Cp(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Cp(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Cp(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::HsEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Hs(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Hs(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Hs(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::LEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().L(idS, YSolid_, p, T);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, rootVSmall);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = 1.0 - YMix[GAS] - YMix[LIQ];

    return massNew;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::radiusForCylinder
(
    const scalar Xi, 
    const scalar V0
)
{
    const scalar a = 2.0*constant::mathematical::pi;
    scalar b = -Xi*constant::mathematical::pi;
    //double c = 0;
    scalar d = -V0;
    scalar root1, root2, root3;
    scalar root;
    
    if (d == 0 || d>0)
    {
        Info<<"ERROR: The volume of cylinder should be positive!! Now the value: "<<V0<<"\n"<<endl;
        return 0;
    } //End if d == 0

    b /= a;
    d /= a;
    
    scalar disc, q, r, dum1, s, t, term1, r13;
    q = (-(b*b))/9.0;
    r = -(27.0*d) + b*(-2.0*(b*b));
    r /= 54.0;
    disc = q*q*q + r*r;
    term1 = (b/3.0);
    
    if (disc > 0) { // one root real, two are complex
        s = r + Foam::sqrt(disc);
        s = ((s < 0) ? -Foam::pow(-s, (1.0/3.0)) : Foam::pow(s, (1.0/3.0)));
        t = r - Foam::sqrt(disc);
        t = ((t < 0) ? -Foam::pow(-t, (1.0/3.0)) : Foam::pow(t, (1.0/3.0)));
        root = -term1 + s + t;
        return root;
    } 
    // End if (disc > 0)
    
    // The remaining options are all real
    if (disc == 0){ // All roots real, at least two are equal.
        r13 = ((r < 0) ? -Foam::pow(-r,(1.0/3.0)) : Foam::pow(r,(1.0/3.0)));
        root1 = -term1 + 2.0*r13;
        root2 = -(r13 + term1);
        root = ((root1 > root2) ? root1 : root2);
        return root;
    } // End if (disc == 0)
    
    // Only option left is that all roots are real and unequal (to get here, q < 0)
    q = -q;
    dum1 = q*q*q;
    dum1 = Foam::acos(r/Foam::sqrt(dum1));
    r13 = 2.0*Foam::sqrt(q);
    root1 = -term1 + r13*Foam::cos(dum1/3.0);
    root2 = -term1 + r13*Foam::cos((dum1 + 2.0*constant::mathematical::pi)/3.0);
    root3 = -term1 + r13*Foam::cos((dum1 + 4.0*constant::mathematical::pi)/3.0);
    if (root1 > root2){
        root = ((root1 > root3) ? root1 : root3);
    }
    else{
        root = ((root2 > root3) ? root2 : root3);
    }
    return root;
}  //End of cubicSolve

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::Area_Sph
(
    const scalar radius
)
{
    return 4.0*constant::mathematical::pi*Foam::sqr(radius);
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::Area_cylinderL
(
    const scalar radius,
    const scalar Xi
)
{
    return 2.0*constant::mathematical::pi*(3.0*Foam::sqr(radius)-radius*Xi);
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::Volume_cylinderL
(
    const scalar radius,
    const scalar Xi
)
{
    return constant::mathematical::pi*(2.0*Foam::pow3(radius)-sqr(radius)*Xi);
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::Vol_Rin 
(
    const scalar rin,
    const scalar rout
)
{
    
    return (4./3.)*constant::mathematical::pi*(Foam::pow3(rout) - Foam::pow3(rin));
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::Vol_Rin_cylinderL
(
    const scalar rin,
    const scalar rout,
    const scalar Xi
)
{
    return constant::mathematical::pi*(2.0*Foam::pow3(rout) - 2.0*Foam::pow3(rin) - (Foam::sqr(rout) - Foam::sqr(rin))*Xi);
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::R_Par
(
    const scalar rin,
    const scalar rout
)
{
    scalar temp ,r;
    
    temp = (Foam::pow3(rin)+Foam::pow3(rout))/2.0;
    
    r = Foam::cbrt(temp);
    
    return r;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::R_Par_cylinderL
(
    const scalar rin,
    const scalar rout,
    const scalar Xi
)
{
    scalar temp ,r;
    
    temp = (2.0*constant::mathematical::pi*(pow3(rout) + pow3(rin))-constant::mathematical::pi*(pow(rout,2.0) + pow(rin,2.0))*Xi)/2.0;
    
    r = radiusForCylinder(Xi, temp);
    
    return r;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::d_dr
(
    const scalar kp,
    const scalar Ab,
    const scalar ri,
    const scalar rj
)
{
    return  kp*Ab*rj/(ri*(ri-rj));  
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::d_dr_cylinderL
(
    const scalar kp,
    const scalar Ab,
    const scalar ri,
    const scalar rj,
    const scalar Xi
)
{
   return  kp*Ab*Xi/((ri*Xi-3.0*Foam::sqr(ri))*Foam::log((3.0-Xi/rj)/(3.0-Xi/ri)));  
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::Fb1 
(
    const scalar Tb1
)
{
        
    return  Foam::pow(10.,(8.07131-1730.63/(Tb1-39.724)))/760.0;
}

template<class ParcelType>
bool Foam::ReactingMultiphaseParcel<ParcelType>::eq4_check
(
    const scalar h_coe,
    const scalar emissi,
    const scalar Ste_Bol,
    const scalar kp3,
    const scalar Ab3,
    const scalar rb3,
    const scalar rp3,
    const scalar Tg,
    const scalar Tp3,
    const scalar G,
    const scalar Source,
    scalar& Tb3
)
{
    //the mothed used here is General formula for roots, see wikipedia, well it is quite unstable in practice.
    scalar aa, dd, ee, delta0, delta, delta1, Q, S, q,sqrt1,sqrt2,sqrt3;
    aa = emissi*Ste_Bol*Ab3;
    if ((rb3-rp3)  == 0)
    {
            FatalErrorInFunction
            <<"(rb3-rp3)  == 0" << exit(FatalError);
    }

    dd = h_coe*Ab3+Ab3*kp3*rp3/(rb3*(rb3-rp3));

    ee = -Ab3*kp3*rp3/(rb3*(rb3-rp3))*Tp3-h_coe*Ab3*Tg-aa*(G/(4.*Ste_Bol))-Source;
//  Info<<"Tg emission:   "<<G/(4.*Ste_Bol)<<endl;

    delta = 256.0*Foam::pow3(aa)*Foam::pow3(ee)-27.0*aa*aa*Foam::pow4(dd);

    sqrt1 = -27.0*delta;
    
// Info<<"debug information ****** checkPoint 1  ******"<<endl;
    if (sqrt1  < 0)
    {
        return false;
    }
    delta1 = 27.0*aa*dd*dd;
    delta0 = 12.0*aa*ee;                
    Q = Foam::cbrt(0.5*(delta1+Foam::sqrt(sqrt1)));
// Info<<"debug information ****** checkPoint 2  ******"<<endl;    
    if (mag(Q) < rootVSmall)
    {
        return false;
    }
    
    sqrt2 = (Q+delta0/Q)/(3.0*aa);

    if (sqrt2  < 0)
    {
        return false;
    }
    
    S = 0.5*Foam::sqrt(sqrt2);
    
    if (mag(S) < small)
    {
        return false;
    } 
    
    if (mag(aa) < rootVSmall)
    {
        return false;
    } 
  
    q = dd/aa;
    
    sqrt3 = (-4*S*S+q/S);
    
    if (sqrt3  < 0)
    {
        return false;
    }
    
    Tb3 = -S+0.5*Foam::sqrt(sqrt3);         
    
    return true;
}
        
//- eq4 Henrik C&F for cylinderL
template<class ParcelType>
bool Foam::ReactingMultiphaseParcel<ParcelType>::eq4_cylinderL_check
(
    const scalar h_coe,
    const scalar emissi,
    const scalar Ste_Bol,
    const scalar kp3,
    const scalar Ab3,
    const scalar rb3,
    const scalar rp3,
    const scalar Tg,
    const scalar Tp3,
    const scalar G,
    const scalar Source,
    const scalar Xi,
    scalar& Tb3
)
{
    scalar aa, dd, ee, delta0, delta, delta1, Q, S, q,sqrt1,sqrt2,sqrt3;
    aa = emissi*Ste_Bol*Ab3;
    dd = h_coe*Ab3+d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi);
    ee = -d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi)*Tp3-h_coe*Ab3*Tg-aa*(G/(4.*Ste_Bol))-Source;
    delta = 256.0*Foam::pow3(aa)*Foam::pow3(ee)-27.0*aa*aa*Foam::pow4(dd);
    sqrt1 = -27.0*delta;
    if (sqrt1  < 0)
    {
        return false;
    }
    delta1 = 27.0*aa*dd*dd;
    delta0 = 12.0*aa*ee;                
    Q = Foam::cbrt(0.5*(delta1+Foam::sqrt(sqrt1)));
    sqrt2 = (Q+delta0/Q)/(3.0*aa);
    if (sqrt2 < 0)
    {
        return false;
    }
    S = 0.5*Foam::sqrt(sqrt2);
    q = dd/aa;
    sqrt3 = (-4*S*S+q/S);
    if (sqrt3 < 0)
    {
        return false;
    }
    
    Tb3 = -S+0.5*Foam::sqrt(-4*S*S+q/S);
    
    return true;
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::eq4 
(
    const scalar h_coe,
    const scalar emissi,
    const scalar Ste_Bol,
    const scalar kp3,
    const scalar Ab3,
    const scalar rb3,
    const scalar rp3,
    const scalar Tg,
    const scalar Tp3,
    const scalar G,
    const scalar Source
)
{
    //the mothed used here is General formula for roots, see wikipedia, well it is quite unstable in practice.
    scalar aa, dd, ee, delta0, delta, delta1, Q, S, q, Tb3,sqrt1,sqrt2,sqrt3;
    aa = emissi*Ste_Bol*Ab3;
    dd = h_coe*Ab3+Ab3*kp3*rp3/(rb3*(rb3-rp3));
    ee = -Ab3*kp3*rp3/(rb3*(rb3-rp3))*Tp3-h_coe*Ab3*Tg-aa*(G/(4.*Ste_Bol))-Source;
//  Info<<"Tg emission:   "<<G/(4.*Ste_Bol)<<endl;
    delta = 256.0*Foam::pow3(aa)*Foam::pow3(ee)-27.0*aa*aa*Foam::pow4(dd);
    delta1 = 27.0*aa*dd*dd;
    delta0 = 12.0*aa*ee;
    sqrt1 = -27.0*delta;
    if (sqrt1  < 0)
    {
            FatalErrorInFunction
            <<"sqrt1 < 0" << exit(FatalError);
    }
                    
    Q = Foam::cbrt(0.5*(delta1+Foam::sqrt(sqrt1)));
    sqrt2 = (Q+delta0/Q)/(3.0*aa);
    if (sqrt2 < 0)
    {
            FatalErrorInFunction
            <<"sqrt2 < 0" << exit(FatalError);
    }
    S = 0.5*Foam::sqrt(sqrt2);
    q = dd/aa;
    sqrt3 = (-4*S*S+q/S);
    if (sqrt3 < 0)
    {
            FatalErrorInFunction
            <<"sqrt3 < 0" << exit(FatalError);
    }               
    Tb3 = -S+0.5*Foam::sqrt(sqrt3);
    
    return Tb3;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::eq4_cylinderL
(
    const scalar h_coe,
    const scalar emissi,
    const scalar Ste_Bol,
    const scalar kp3,
    const scalar Ab3,
    const scalar rb3,
    const scalar rp3,
    const scalar Tg,
    const scalar Tp3,
    const scalar G,
    const scalar Source,
    const scalar Xi
)
{
    scalar aa, dd, ee, delta0, delta, delta1, Q, S, q, Tb3;
    aa = emissi*Ste_Bol*Ab3;
    dd = h_coe*Ab3+d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi);
    ee = -d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi)*Tp3-h_coe*Ab3*Tg-aa*(G/(4.*Ste_Bol))-Source;
    delta = 256.0*Foam::pow3(aa)*Foam::pow3(ee)-27.0*aa*aa*Foam::pow4(dd);
    delta1 = 27.0*aa*dd*dd;
    delta0 = 12.0*aa*ee;
    Q = Foam::cbrt(0.5*(delta1+Foam::sqrt(-27.0*delta)));
    S = 0.5*Foam::sqrt((Q+delta0/Q)/(3.0*aa));
    q = dd/aa;
    
    Tb3 = -S+0.5*Foam::sqrt(-4*S*S+q/S);
    
    return Tb3;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::eq4_Explicit
(
    const scalar h_coe,
    const scalar emissi,
    const scalar Ste_Bol,
    const scalar kp3,
    const scalar Ab3,
    const scalar rb3,
    const scalar rp3,
    const scalar Tg,
    const scalar Tp3,
    const scalar Tb3_old,
    const scalar G,
    const scalar Source
)
{
    
// Info<<"debug information ****** eq4_Explicit called ******"<<endl;
    //the mothed used here is General formula for roots, see wikipedia, well it is quite unstable in practice.
    scalar aa, dd, ee, Tb3;
    aa = emissi*Ste_Bol*Ab3;
    dd = h_coe*Ab3+Ab3*kp3*rp3/(rb3*(rb3-rp3));
    ee = -Ab3*kp3*rp3/(rb3*(rb3-rp3))*Tp3+aa*Foam::pow4(Tb3_old)-h_coe*Ab3*Tg-aa*(G/(4.*Ste_Bol))-Source;

// Info<<"debug information ****** Tp3: "<<Tp3<< "******"<<"Tg: "<<Tg<< "******"<<endl;
// Info<<"ee1: "<<-Ab3*kp3*rp3/(rb3*(rb3-rp3))*Tp3<< "******"<<"ee2: "<<aa*Foam::pow4(Tb3_old)<< "******"<<"ee3: "<<-h_coe*Ab3*Tg<< "******"<<"ee4: "<<-aa*(G/(4.*Ste_Bol))-Source<< "******"<<"ee5: "<<-Source<< "******"<<endl;
//     
    Tb3 = -ee/dd;
// Info<<"debug information ****** eq4_Explicit result ******"<<Tb3<<endl;   
    return Tb3;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::eq4_cylinderL_Explicit
(
    const scalar h_coe,
    const scalar emissi,
    const scalar Ste_Bol,
    const scalar kp3,
    const scalar Ab3,
    const scalar rb3,
    const scalar rp3,
    const scalar Tg,
    const scalar Tp3,
    const scalar Tb3_old,
    const scalar G,
    const scalar Source,
    const scalar Xi
)
{
    scalar aa, dd, ee, Tb3;
    aa = emissi*Ste_Bol*Ab3;
    dd = h_coe*Ab3+d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi);
    ee = -d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi)*Tp3+aa*Foam::pow4(Tb3_old)-h_coe*Ab3*Tg-aa*(G/(4.*Ste_Bol))-Source;
    
    Tb3 = -ee/dd;
    
//     Info<<"convection: "<<h_coe*Ab3*(Tb3 - Tg)<<endl;
//     Info<<"conduction: "<<d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi)*(Tb3 - Tp3)<<endl;
//     Info<<"radiation: "<<aa*(Foam::pow4(Tb3_old) - G/(4.*Ste_Bol))<<endl;
//     Info<<"Source: "<<Source<<endl;
//     Info<<"Tb3: "<<Tb3<<endl;
//     Info<<"Tp3: "<<Tp3<<endl;
//     
    return Tb3;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::eq7_2
(
    const scalar kp1,
    const scalar kp2,
    const scalar Ab1,
    const scalar rb1,
    const scalar rp1,
    const scalar rp2,
    const scalar Tp1,
    const scalar Tp2,
    const scalar Qb1    
)
{
    scalar temp1 = d_dr(kp2,Ab1,rb1,rp2), temp2 = d_dr(kp1,Ab1,rb1,rp1);
        
    return (-Qb1+temp1*Tp2-temp2*Tp1)/(temp1-temp2);    
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::eq7_2_cylinderL
(
    const scalar kp1,
    const scalar kp2,
    const scalar Ab1,
    const scalar rb1,
    const scalar rp1,
    const scalar rp2,
    const scalar Tp1,
    const scalar Tp2,
    const scalar Qb1,
    const scalar Xi
)
{
    scalar temp1 = d_dr_cylinderL(kp2,Ab1,rb1,rp2,Xi), temp2 = d_dr_cylinderL(kp1,Ab1,rb1,rp1,Xi);
    
// Info<<"source: "<<Qb1<<endl;
// Info<<"(temp1-temp2): "<<(temp1-temp2)<<endl;
// Info<<"source contribution: "<<Qb1/(temp1-temp2)<<endl;

    return (-Qb1+temp1*Tp2-temp2*Tp1)/(temp1-temp2);    
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::eq7_3
(
    const scalar kp1,
    const scalar kp2,
    const scalar Ab1,
    const scalar rb1,
    const scalar rp1,
    const scalar rp2,
    const scalar Tp1,
    const scalar Tp2,
    const scalar Fb    
)
{
    scalar temp1 = d_dr(kp2,Ab1,rb1,rp2), temp2 = d_dr(kp1,Ab1,rb1,rp1);
        
    return (temp1*(1.-Fb)*Tp2-temp2*Tp1)/(temp1*(1.-Fb)-temp2);   
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::eq7_3_cylinderL
(
    const scalar kp1,
    const scalar kp2,
    const scalar Ab1,
    const scalar rb1,
    const scalar rp1,
    const scalar rp2,
    const scalar Tp1,
    const scalar Tp2,
    const scalar Fb,
    const scalar Xi
)
{
    scalar temp1 = d_dr_cylinderL(kp2,Ab1,rb1,rp2,Xi), temp2 = d_dr_cylinderL(kp1,Ab1,rb1,rp1,Xi);
        
    return (temp1*(1.-Fb)*Tp2-temp2*Tp1)/(temp1*(1.-Fb)-temp2);   
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::cp_p 
(
    const label layer,
    const scalar Tp,
    const scalar moist_WB
)
{ 
    scalar heat_capacity=0.;
    
    scalar Tp_2 = Tp;
//     if (Tp_2>=400)
//     {
//         Tp_2=400;
//     }
//     
    if (layer == 0)
    {
        /* moist wood */
        scalar c_p_wood = 4.206*Tp_2 - 37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
        scalar A = 1.0e3*((0.02355*Tp_2 - 1.320*moist_WB/(1.0 - moist_WB) - 6.191)*moist_WB/(1.0 - moist_WB));
        heat_capacity = c_p_wood*(1.0 - moist_WB) + 4185.5*moist_WB + A;    
    }
    else if (layer == 1)
    {
        /* dry wood */
        heat_capacity = 4.206*Tp_2 - 37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
    }
    else if (layer == 2)
    {
        /* char */
        heat_capacity = -334.0 + 4410.0e-3*Tp_2 - 3160.0e-6*Foam::pow(Tp_2,2.0) + 1010.0e-9*Foam::pow(Tp_2,3.0) - 119.0e-12*Foam::pow(Tp_2,4.0); /* Table 5 Thunman et al., Energy & Fuels 2001 */
    }
    else if (layer == 3)
    {
        heat_capacity = 754.0 + 0.586*(Tp_2-273.15);     /* Sven */
    }
    
    return heat_capacity;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::cp_p_modified
(
    const label layer,
    const scalar Tp,
    const scalar moist_WB
)
{ 
    scalar heat_capacity=0.;
    
    scalar Tp_2 = Tp;
    
    if (Tp_2 < 273)
    {
        Tp_2=273;
    }
    
    if (layer == 0)
    {
        /* moist wood */
        scalar c_p_wood = 4.206*Tp_2 - 37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
        scalar A = 1.0e3*((0.02355*Tp_2 - 1.320*moist_WB/(1.0 - moist_WB) - 6.191)*moist_WB/(1.0 - moist_WB));
        heat_capacity = c_p_wood*(1.0 - moist_WB) + 4185.5*moist_WB + A;    
    }
    else if (layer == 1)
    {
        /* dry wood */
        heat_capacity = 4.206*Tp_2 - 37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
    }
    else if (layer == 2)
    {
        /* char */
        heat_capacity = 420+2.09*Tp_2+0.000685*Foam::pow(Tp_2,2.0); /* Ramin mehrabian., Fuels process techology 2012 */
    }
    else if (layer == 3)
    {
        heat_capacity = 754.0 + 0.586*(Tp_2-273.15);     /* Sven */
    }
    
    return heat_capacity;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::rho_p 
(
    const label layer,
    const scalar Tp,
    const scalar moist_WB
)
{
    scalar density=0;
    
    if (layer == 0)
    {
        density = 750.0;
    }
    else if (layer == 1)
    {
        density = (1.0 - moist_WB)*750.0; 
    }
    else if (layer == 2)
    {
        density = 150.0;
    }
    else if (layer == 3)
    {
        scalar dgas = 101325.0*(2.0*14.00672*1.0e-3)/(8.3145*Tp);
        density = 2000.0*(1.0-0.65) + 0.65*dgas;          /* Sven */
    }
    
    return density;
}
    
template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::kp_p 
(
    const label layer,
    const scalar Tp
)
{
    scalar heat_conductivity=0;
    
    if (layer == 0)
    {
        heat_conductivity = 0.2;                 
    }
    else if (layer == 1)
    {
        heat_conductivity = 0.11;
    }
    else if (layer == 2)
    {
        heat_conductivity = 0.052;
    }
    else if (layer == 3)
    {
        scalar kgas = (-0.000000000003493)*Foam::pow3(Tp) + 0.000000003319264*Foam::sqr(Tp) + 0.000060059499759*Tp + 0.008533051948052;
        heat_conductivity = 1.03*(1.0-0.65) + 0.65*kgas;  /* Sven */
    }
    
    return heat_conductivity;
} 

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::kpp 
(
    const scalar Tp
)
{
    scalar heat_conductivity=0;
    
    label layer = 3;
    
    if (layer == 0)
    {
        heat_conductivity = 0.2;                 
    }
    else if (layer == 1)
    {
        heat_conductivity = 0.11;
    }
    else if (layer == 2)
    {
        heat_conductivity = 0.052;
    }
    else if (layer == 3)
    {
        scalar kgas = (-0.000000000003493)*Foam::pow3(Tp) + 0.000000003319264*Foam::sqr(Tp) + 0.000060059499759*Tp + 0.008533051948052;
        heat_conductivity = 1.03*(1.0-0.65) + 0.65*kgas;  /* Sven */
    }
    
    return heat_conductivity;
}    
    

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::cp_water_vapor 
(
    const scalar T
)
{
    scalar A, B, C, D, E, Temp;
    A = 30.09200;
    B = 6.832514;
    C = 6.793435;
    D = -2.534480;
    E = 0.082139;
    Temp = T/1000.;
    return (A+B*Temp+C*Foam::sqr(Temp)+D*Foam::pow3(Temp)+E/Foam::sqr(Temp))/18.*1000.;    
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::deltaHvap 
(
    const scalar T
)
{
    scalar A, B, C, dH;
    A = 51798.3515625;
    B = -10.6430501938;
    C = -0.0515099391341;
    if (T < 273.15)
    {
        dH = 45.054e3;
    }
    else if ( T > 433.15 )
    {
        dH = 37.518e3;
    }
    else
    {
        dH = A + B*T + C*Foam::sqr(T);
    }
    
    return dH/18.0e-3; /* convert J/mol to J/kg */
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Re-use correction from reacting parcel
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();
        
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    canCombust_ = 1;
    const label particleShape = cloud.constProps().parcelShape();
    scalar Xi;
    
    if (particleShape == 1)
    {
        Xi = cloud.constProps().xi0();
    }
    else
    {
        Xi = 1;
    }
    scalar equivalent_d;
    QLayer_ = 0.0;
    
    const bool coarseGrid = cloud.constProps().coarseGrid();
    const bool coarseGridSurfaceCombustion = cloud.constProps().coarseGridSurfaceCombustion();
    const bool particleToGasOneWayHtc = cloud.constProps().disableHtc();
    const bool noCombustion = cloud.constProps().noCombustion();
    const bool heatRatio = cloud.constProps().heatRatioFlag();
    const bool usePpConduction = cloud.constProps().ppConductionFlag();
    
    const bool limitedCombustionFlag = cloud.constProps().controledCombustion();
    const scalar combustionCriterion = cloud.constProps().combustionAt();
    const label shrinkageMode = cloud.constProps().shrinkageMode();
    bool noCombustionDueToDevo = false;
    
    const bool surfaceExplicitFlag = cloud.constProps().surfaceExplicit();
    
    const scalar radiationCorrection = cloud.constProps().radiationCorrection();
    const bool usePpRadiation = cloud.constProps().ppRadiationFlag();
    
    const scalar np0 = this->nParticle_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
//     const scalar mass0 = this->mass();

    const scalar pc = td.pc();
    
    // mass at beginning [kg]
    scalar massit = mp0_+mp1_+mp2_+mp3_;
    const scalarField& YMix = this->Y_;
//     const label idG = composition.idGas();
    const label idL = composition.idLiquid();
    const label idS = composition.idSolid();
    const label idDryWood= composition.localId(idS,"wood");
    const label idActiveDryWood= composition.localId(idS,"activeDryWood");
    const label idChar= composition.localId(idS,"C");
    const label idAsh= composition.localId(idS,"ash");
    const label idwater= composition.localId(idL,"H2O");
    
    // initial fraction
    const scalar Ywood00 = composition.YMixture0()[idS]*composition.Y0(idS)[idDryWood];
    const scalar Yash00 = composition.YMixture0()[idS]*composition.Y0(idS)[idAsh];
    const scalar Ywater00 = composition.YMixture0()[idL];
    const scalar YashDB00 = Yash00/(1.-Ywater00);
    const scalar YdryDB00 = 1.-Yash00/(1.-Ywater00);
    const scalar ratioWoodMoist = Ywood00/Ywater00;
//    const scalar ashDryWood = (Yash00 + Ywood00)/Ywood00;

    // ash content in char layer
    scalar ash_inchar = ash_inchar_t_;
    scalar ash_inchar_old = ash_inchar_t_;    

// SC * * * * * * * * * * * * * * * srinkage factor need to be adjusted according to the particle type  * * * * * * * * * * * * * * * * * * * * * * // 
    // srinkage factor TODO make it read from const dict
    const scalar drySrin = 0.10;
    //const scalar devoSrin = 0.28;
    const scalar devoSrin = 0.39;
    scalar charSrinHardCoded = 0.95;
    if (shrinkageMode == 0)
    {
        charSrinHardCoded = cloud.constProps().shrinkageFactorAsh();
    }
    else
    {
        charSrinHardCoded = 0.95;
    }
    const scalar charSrin = charSrinHardCoded;
// -------------------------srinkage factor need to be adjusted according to the particle type---

// limit of the mass left
    const scalar p_WoodLeft = 0.003333;
    
    // just for loop
    label i;

    //force update fraction
    this->Y_[GAS] = 0.0;
    this->Y_[LIQ] = mp0_*Ywater00/massit;
    this->Y_[SLD] = 1.0 - this->Y_[GAS] - this->Y_[LIQ];
 

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);
    scalarField dMassPCSoild(YSolid_.size(), 0.0);
    
    // Mass transfer due to devolatilisation
    scalarField dMassDV(YGas_.size(), 0.0);
    scalarField dMassDVTemp(YGas_.size(), 0.0);
    
    // Solid mass changes due to devolatilisation
    scalarField dMassDVSoild(YSolid_.size(), 0.0);
    scalarField dMassDVSoildTemp(YSolid_.size(), 0.0);

    // Change in carrier phase composition due to surface reactions
    scalarField dMassSRGas(YGas_.size(), 0.0);
    
    //readlly needed?? TODO
    scalarField dMassSRGasvid_(7, 0.0);//debug veriable
    
    
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField dMassSRCarrier(composition.carrier().species().size(), 0.0);
    scalarField dHeatSRCarrier(8, 0.0);//size according to surfacereaction model, here have four model

    // Change in carrier phase composition due to surface reactions in total
    scalarField dMassSRCarrierTot(composition.carrier().species().size(), 0.0);
    
    // Zero field for mass transfer
    scalarField dMassGasZero(YGas_.size(), 0.0);
    scalarField dMassLiquidZero(YLiquid_.size(), 0.0);
    
    // Mass transfer in gas, liquid, and solid phase
    scalarField dMassGas(YGas_.size(), 0.0);
    scalarField dMassLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSolid(YSolid_.size(), 0.0);

    scalar Fb1Value = 0;

    // Calc heat transfer coefficient [W K-1]
//      scalar h_coe=50, emissi=0.8;
    scalar h_coe;
    // constant
    // Ste_Bol      W m-2 K-4   Stefan–Boltzmann constant
    // emissi       -       emissivity
    // Lat_Heat_wat J kg-1  Latent Heat Water
    // R_gas        J K−1 mol−1 gas constant 
    const scalar emissi = cloud.constProps().epsilon0(), Ste_Bol = physicoChemical::sigma.value(), /*Lat_Heat_wat = 2.26E06,*/ T_boiling = 373.15 ; 
    /*h_coe = cloud.constProps().hCoeff(), */
    
    // solver related parameters
    // t_step time step [s], t_end conversion time [s]
    // tolerance tolerance level for tempearture [-], t_current current time [s]
    // tolerance_Tb tolerance level during the calculation [-]
    // tolerance_Tp tolerance level during the calculation [-]
    // maxIters maximum numbers of interation for Tp [-]  
    scalarField tolerance_Tp(4), tolerance_Tb(4);

    scalar t_current = 0, tolerance = 1.0e-07;
    label maxIters = 100;   
     
    scalar t_step = cloud.constProps().deltaTime();
    
    // rp particle layer diameter [m]
    // Ab particle boundary surface area [m2]
    // Vp layer volume [m3]
    // mp layer mass [kg]
    // dTp/dt [K/s]
    // some intermidiete value 
    scalarField Vp(4), mp(4), mp_new(4), mp_old(4), Vp_old(4), Ab(4), dTpdt(4), cp(4), rb(4), rhop(4),
                rp(4), alphap(4),Tp_old(4), Tp_lastiter(4),Tb_lastiter(4), Tp(4), Tb(4), Tb_old(4), kp(4);


      //- Initialise layer variables at start, only once!!
    if (this->age_ == 0.)
    {
        
        massit = mp0_+mp1_+mp2_+mp3_;     
         //force update fraction
        this->Y_[GAS] = 0.0;
        this->Y_[LIQ] = mp0_*Ywater00/massit;
        this->Y_[SLD] = 1.0 - this->Y_[GAS] - this->Y_[LIQ];       
        //Force correcting solid mass fraction
        this->mass0_ = massit;         
    }
                
    Tp[0] = this->Tp0_;
    Tp[1] = this->Tp1_;
    Tp[2] = this->Tp2_;
    Tp[3] = this->Tp3_;
    
    Tb[0] = this->Tb0_;
    Tb[1] = this->Tb1_;
    Tb[2] = this->Tb2_;
    Tb[3] = this->Tb3_;

    rb[0] = this->rb0_;
    rb[1] = this->rb1_;
    rb[2] = this->rb2_;
    rb[3] = this->rb3_;

    if (particleShape == 1)
    {
        rp[0] = R_Par_cylinderL(0.,rb[0],Xi);
        rp[1] = R_Par_cylinderL(rb[0],rb[1],Xi);
        rp[2] = R_Par_cylinderL(rb[1],rb[2],Xi);
        rp[3] = R_Par_cylinderL(rb[2],rb[3],Xi); 
        
        Ab[0] = Area_cylinderL(rb[0],Xi);
        Ab[1] = Area_cylinderL(rb[1],Xi);
        Ab[2] = Area_cylinderL(rb[2],Xi);
        Ab[3] = Area_cylinderL(rb[3],Xi);

        Vp[0] = Volume_cylinderL(rb[0],Xi);
        Vp[1] = Vol_Rin_cylinderL(rb[0],rb[1],Xi);
        Vp[2] = Vol_Rin_cylinderL(rb[1],rb[2],Xi);
        Vp[3] = Vol_Rin_cylinderL(rb[2],rb[3],Xi);
    }
    else
    {
        rp[0] = R_Par(0.,rb[0]);
        rp[1] = R_Par(rb[0],rb[1]);
        rp[2] = R_Par(rb[1],rb[2]);
        rp[3] = R_Par(rb[2],rb[3]); 
        
        Ab[0] = Area_Sph(rb[0]);
        Ab[1] = Area_Sph(rb[1]);
        Ab[2] = Area_Sph(rb[2]);
        Ab[3] = Area_Sph(rb[3]);

        Vp[0] = (4.0/3.0)*constant::mathematical::pi*Foam::pow3(rb[0]);
        Vp[1] = (4.0/3.0)*constant::mathematical::pi*(Foam::pow3(rb[1])-Foam::pow3(rb[0]));
        Vp[2] = (4.0/3.0)*constant::mathematical::pi*(Foam::pow3(rb[2])-Foam::pow3(rb[1]));
        Vp[3] = (4.0/3.0)*constant::mathematical::pi*(Foam::pow3(rb[3])-Foam::pow3(rb[2]));
    }
 
    mp[0] = this->mp0_;
    mp[1] = this->mp1_;
    mp[2] = this->mp2_;
    mp[3] = this->mp3_;

    tolerance_Tb[0] = GREAT;
    tolerance_Tb[1] = GREAT;
    tolerance_Tb[2] = GREAT;
    tolerance_Tb[3] = GREAT;
    tolerance_Tp[0] = GREAT;
    tolerance_Tp[1] = GREAT;
    tolerance_Tp[2] = GREAT;
    tolerance_Tp[3] = GREAT;
    
    // Some variables used to calculate heat transfer between layers 
    scalar Qp0to1, Qp1to2, Qp2to3, Qg0to1, Qg1to2, /*Qg2to3,*/ mp0to1, mp1to2, mp2to3, mGas0, mGas1, cpG00to1_temp, 
           cpG01to2_temp, /*cpG02to3_temp, */cpG11to2_temp, /*cpG12to3_temp,*/ cpP0to1_temp, cpP1to2_temp, cpP2to3_temp;


    //Force correcting solid mass fraction
    YSolid_[idDryWood] = (mp[0])*Ywood00/(YMix[SLD]*massit);
    YSolid_[idActiveDryWood] = mp[1]*composition.Y0(idS)[idDryWood]/(YMix[SLD]*massit);
    YSolid_[idAsh] = this->mass0_*Yash00/(YMix[SLD]*massit);
    YSolid_[idChar] =  1.0 - YSolid_[idDryWood] - YSolid_[idActiveDryWood] - YSolid_[idAsh];
           
 
    //Henrik layer for devolatilisation
    scalar numberSubLayer = 5.;//TODO this value should read from the properties
    label subLayer;
    scalar deltRadiusSubLayer;
    scalar deltTemperature1;
    scalar b0Temperature1;
    scalar deltTemperature2;
    scalar b0Temperature2;    
    scalar radiusSubLayerL;
    scalar radiusSubLayerM;
    scalar radiusSubLayerH;
    scalar massSubLayer;
    scalar temperatureSubLayer;
    
    scalar QDevo = 0.0;

    //Radiation
    tetIndices tetIs = this->currentTetIndices();
    
    scalar GcInitial = Foam::pow4(273.0)*(4.*Ste_Bol);
    
    if (cloud.radiation())
    {
        GcInitial = td.GInterp().interpolate(this->coordinates(), tetIs);
        if (GcInitial < Foam::pow4(200.0)*(4.*Ste_Bol))
        {
            GcInitial = Foam::pow4(273.0)*(4.*Ste_Bol);
        }
    }
    
    
//     const scalar Gc = 601275.983376194, TInifinit =1050.; //Lu
//     const scalar Gc = 329670.829, TInifinit =1098.;//Wurzenberger
//    const scalar TInifinit = this->T_;
    scalar TInifinit = td.Tc();//td.TInterp().interpolate(this->coordinates(), tetIs);
            // in Thermoparcel
            //- Return const access to the interpolator for continuous
            //  phase temperature field
            //inline const interpolation<scalar>& TInterp() const;
    
    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas, Cpcs;
    
    scalar Res;
    
    // bed voidage
    scalar e_bed;
    //minimum bed voidage
    const scalar e_bedMin = (1.0 - cloud.constProps().alphaMax());

    word UScheme = cloud.solution().interpolationSchemes().lookup("U");
    word TScheme = cloud.solution().interpolationSchemes().lookup("T");
    word GScheme = cloud.solution().interpolationSchemes().lookup("G");

    // Sources
    //~~~~~~~~

//     // Explicit momentum source for particle
//     vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;

//     VID1_ = 0;
//     VID2_ = 0;
    VID3_ = 0;
    VID4_ = 0;
    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Sum Ni*Cpi*Wi of emission species
//     scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(composition.carrier().species().size(), 0.0);
    
     // Average carriers surface concentrations and bed voidage
    // 0-O2, 1-h2o, 2-co2, 3-h2, 4-e_bed
    scalarField Csavg(6, 0.0); 
    Csavg[0] = cloud.o2avg()[this->cell()];
    Csavg[1] = cloud.h2oavg()[this->cell()];
    Csavg[2] = cloud.co2avg()[this->cell()];
    Csavg[3] = cloud.h2avg()[this->cell()];
// Info<<"new particle"<<endl;
    //particle submodel
    do
    {
// Info<<"particle new time step"<<endl;
        //save data last time step
        for ( i=0;i<=3;i++)
        {
            Tp_old[i] = Tp[i];
            Tb_old[i] = Tb[i];
            Tp_lastiter[i] = Tp[i];
            Tb_lastiter[i] = Tb[i];
            mp_old[i] = mp[i];
            mp_new[i] =mp[i];
            Vp_old[i] = Vp[i];            
        }
        ash_inchar_old = ash_inchar;
        //iteration number 
        label iter_0 = 0;
        scalar Vp1_temp;
        scalar Vp2_temp;
        
        //last iteration flag
        bool lastIteration = false;
            
        //iteration loop

        if (coarseGrid)
        {
            // Calc surface values
            if (TScheme == "cellPoint")
            {
                TInifinit = cloud.T()[this->cell()];
            }
            else
            {
                TInifinit = cloud.Tavg()[this->cell()];
            }

            // Surface temperature using two thirds rule
            Ts = (2.0*T0 + TInifinit)/3.0;

            if (Ts < cloud.constProps().TMin())
            {
                if (debug)
                {
                    WarningInFunction
                        << "Limiting parcel surface temperature to "
                        << cloud.constProps().TMin() <<  nl << endl;
                }

                Ts = cloud.constProps().TMin();
            }

            // Assuming thermo props vary linearly with T for small d(T)
            const scalar TRatio = TInifinit/Ts;

            rhos = cloud.rhoavg()[this->cell()]*TRatio;

            Cpcs = cloud.cpavg()[this->cell()];
            mus = cloud.muavg()[this->cell()];
            kappas = cloud.kappaavg()[this->cell()];

            Prs = Cpcs*mus/kappas;
            Prs = max(rootVSmall, Prs);
            
            // Reynolds number
            if (UScheme == "cellPoint")
            {
                Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
            }
            else
            {
                Res = this->Re(rhos, U0, cloud.Uavg()[this->cell()], this->d_, mus);
            }
            
            e_bed = max(cloud.alphaavg()[this->cell()], e_bedMin);
            
             if (cloud.radiation())
            {
                if (GScheme == "cell")
                {
                    GcInitial = cloud.mesh().objectRegistry::template
                        lookupObject<volScalarField>("Gavg")[this->cell()];
                }
                else
                {
                    GcInitial = td.GInterp().interpolate(this->coordinates(), tetIs);
                }
                
                if (GcInitial < Foam::pow4(200.0)*(4.*Ste_Bol))
                {
                    GcInitial = Foam::pow4(T0)*(4.*Ste_Bol);
                }
            }
            
        }
        else
        {
            // Calc surface values
            this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
            
            // Reynolds number
            Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
            
            e_bed = max(1 - (Vp[0]+Vp[1]+Vp[2]+Vp[3])/cloud.pMesh().cellVolumes()[this->cell()], e_bedMin);
        }
        
        //set a limit of maximum particle environment temperature, because at wall region, the gas phase temperature can be too high due to source terms(maybe?) !!!!!!!!!!!! need to be fixed TODO
        scalar maxGcInitial = 4.*Ste_Bol*pow4(1550.0);
        if (GcInitial > maxGcInitial)
        {
            GcInitial = maxGcInitial;
        }
        
         if (TInifinit > 1600.0)
        {
            TInifinit = 1600.0;
        }

//         Info<<"dBeforeDem: "<<this->dBeforeDEM_<<endl;
        const scalar contactedNeighbourNum = this->neighbourNum_; 
        scalar ppGc = GcInitial;
        if(usePpRadiation && contactedNeighbourNum > 0.99)
        {
            if (radiationCorrection < 0)
            {
                ppGc = (1.0 - e_bed)*4*Ste_Bol*this->ppT4sum_/contactedNeighbourNum + e_bed*4.*Ste_Bol*pow4(TInifinit);
            }
            else
            {
                ppGc = 4*Ste_Bol*this->ppT4sum_/contactedNeighbourNum;
            }
            
        
        }
        
        
//         const scalar Gc = GcInitial;
        
        Csavg[4] = e_bed;
        
        if (shrinkageMode == 0)
        {
            //ash layer porosity, ash true density is 2000 kg/m3
            Csavg[5] = 1.0 - mp[3]/Vp[3]/2000.0;
            
//             Info<<"shrinkageMode == 0 here, porosity is : "<<Csavg[5]<<endl;
        }
        else if(shrinkageMode == 1)
        {
//             Info<<"shrinkageMode == 1 here"<<endl;
            //assume porosity is constant in ash is 0.65
            Csavg[5] = 0.65;
        }
        else
        {
            //assume porosity is constant in ash is 0.65
            Csavg[5] = 0.65;
        }

        // Sensible enthalpy transfer from the particle to the carrier phase during  a particle time step   
        scalar dhsTrans_dt = 0.0;    
        
        scalar QCombForLayer = 0.0; //this a a nurmerical stable method, if Qcomb too large, the temperature will be abnormal, instead of add Qcomb to boundary, add to p2 layer.
        
        // Molar flux of species emitted from the particle (kmol/m^2/s)
        scalar Ne = 0.0;
        
        Ne = (rDevo_ - rChar_)/(Ab[3]*t_step*55.0) + rDry_/(Ab[3]*t_step*18.0);
        
        scalar aa = emissi*Ste_Bol*Ab[3];
       
        //iteration loop
        do
        {
            
            // Molar flux of species emitted from the particle (kmol/m^2/s), temp variables duing iteration
            scalar Ne_temp = 0.0;

            //update mass at t+dt
            mp_new[0] = mp_old[0]-rDry_/Ywater00;
            mp_new[1] = mp_old[1]+rDry_*ratioWoodMoist-rDevo_/YdryDB00;
            mp_new[2] = mp_old[2]+rDevo_/YdryDB00*YashDB00+rChar_-rComb_;
            mp_new[3] = mp_old[3]+(ash_inchar_old/mp_old[2])*rComb_;  
            
            if (mp_new[0] > this->mass0_*p_WoodLeft*1.5)
            {
                flagBoiling_ = 1;
            }
            
             if (mp_new[1] > this->mass0_*p_WoodLeft*1.5)
            {
                flagDevo_ = 1;
            }
            
            //update volume, radius, area at T+dt
            Vp[0] = (Vp_old[0]/mp_old[0])*mp_new[0];
            Vp1_temp = (Vp_old[1]/mp_old[1])*(mp_old[1]-rDevo_/YdryDB00);
            Vp[1] = Vp1_temp+(1. - drySrin)*(Vp_old[0] - Vp[0]);
            Vp2_temp = (Vp_old[2]/mp_old[2])*(mp_old[2]-rComb_);
            Vp[2] = Vp2_temp + (1. - devoSrin)*(Vp_old[1] - Vp1_temp);
            Vp[3] = Vp_old[3] + (1. - charSrin)*(Vp_old[2] - Vp2_temp);
// Info<<"start iteration--"<<endl;
// Info<<"mp_old[0]: "<< mp_old[0]<<" mp_new[0]: "<< mp_new[0]<<endl; 
// Info<<"mp_old[1]: "<< mp_old[1]<<" mp_new[1]: "<< mp_new[1]<<endl;
// Info<<"mp_old[2]: "<< mp_old[2]<<" mp_new[2]: "<< mp_new[2]<<endl; 
// Info<<"mp_old[3]: "<< mp_old[3]<<" mp_new[3]: "<< mp_new[3]<<endl; 
// Info<<"Vp[0]"<< Vp[0]<<endl; 
// Info<<"Vp[1]"<< Vp[1]<<endl; 
// Info<<"Vp[2]"<< Vp[2]<<endl; 
// Info<<"Vp[3]"<< Vp[3]<<endl;

            if (particleShape == 1)
            {
                rb[0] = radiusForCylinder(Xi,Vp[0]);
                rp[0] = R_Par_cylinderL(0.0,rb[0],Xi);
                Ab[0] = Area_cylinderL(rb[0],Xi);
                rb[1] = radiusForCylinder(Xi,Vp[0]+Vp[1]);;
                rp[1] = R_Par_cylinderL(rb[0],rb[1],Xi);
                Ab[1] = Area_cylinderL(rb[1],Xi);
                rb[2] = radiusForCylinder(Xi,Vp[0]+Vp[1]+Vp[2]);;
                rp[2] = R_Par_cylinderL(rb[1],rb[2],Xi);
                Ab[2] = Area_cylinderL(rb[2],Xi);
                rb[3] = radiusForCylinder(Xi,Vp[0]+Vp[1]+Vp[2]+Vp[3]);;
                rp[3] = R_Par_cylinderL(rb[2],rb[3],Xi);
                Ab[3] = Area_cylinderL(rb[3],Xi);    
            }
            else
            {
                rb[0] = Foam::cbrt(0.75*Vp[0]/constant::mathematical::pi);
                rp[0] = R_Par(0.0,rb[0]);
                Ab[0] = Area_Sph(rb[0]);
                rb[1] = Foam::cbrt(0.75*(Vp[0]+Vp[1])/constant::mathematical::pi);
                rp[1] = R_Par(rb[0],rb[1]);
                Ab[1] = Area_Sph(rb[1]);
                rb[2] = Foam::cbrt(0.75*(Vp[0]+Vp[1]+Vp[2])/constant::mathematical::pi);
                rp[2] = R_Par(rb[1],rb[2]);
                Ab[2] = Area_Sph(rb[2]);
                rb[3] = Foam::cbrt(0.75*(Vp[0]+Vp[1]+Vp[2]+Vp[3])/constant::mathematical::pi);
                rp[3] = R_Par(rb[2],rb[3]);
                Ab[3] = Area_Sph(rb[3]);    
            }
            
            //update properties
            rhop[0] = mp_new[0]/Vp[0];//rho_p(0,Tp[0],Ywater00);
            rhop[1] = mp_new[1]/Vp[1];//rho_p(1,Tp[1],Ywater00);
            rhop[2] = mp_new[2]/Vp[2];//rho_p(2,Tp[2],Ywater00);
            rhop[3] = mp_new[3]/Vp[3];//rho_p(3,Tp[3],Ywater00);

            kp[0] = kp_p(0,Tp[0]);
            kp[1] = kp_p(1,Tp[1]);
            kp[2] = kp_p(2,Tp[2]);
            kp[3] = kp_p(3,Tp[3]);
            
            cp[0] = cp_p_modified(0,Tp[0],Ywater00);
            cp[1] = cp_p_modified(1,Tp[1],Ywater00);
            cp[2] = cp_p_modified(2,Tp[2],Ywater00);
            cp[3] = cp_p_modified(3,Tp[3],Ywater00);

            alphap[0] = kp[0]/(rhop[0]*cp[0]);
            alphap[1] = kp[1]/(rhop[1]*cp[1]);
            alphap[2] = kp[2]/(rhop[2]*cp[2]);
            alphap[3] = kp[3]/(rhop[3]*cp[3]);  
            
           if (alphap[0] < 0)
           {
               alphap[0] = 1.0e-7;
           }
            
           if (alphap[1] < 0)
           {
               alphap[1] = 7.0e-8;
           }
           
            if (alphap[2] < 0)
           {
               alphap[2] = 3.0e-7;
           }
            // storage Tb at last iteration  
            for ( i=0;i<=3; i++)
            {
                Tb_lastiter[i] = Tb[i];         
            }            
            
            //Heat transfer coefficient
            //h_coe = cloud.heatTransfer().htc(2.0*rb[3], Res, Prs, kappas, NCpW);
            //TODO need proper correlations
            //h_coe = kappas/(this->d_)*(2.0+1.1*cbrt(Prs)*pow(Res,0.6));
            h_coe = kappas/(this->d_)*(1.77+0.29*pow(e_bed,-0.81)*sqrt(Prs)*pow(Res,0.73));

            scalar correctedEmissi = emissi;
            if (flagBoiling_ == 1 || flagDevo_ == 1)
            {
                if(usePpRadiation && contactedNeighbourNum > 0.99)
                {
                    GcInitial = ppGc;
                }
                else
                {
                    if (radiationCorrection < 0)
                    {
                        correctedEmissi = emissi;
                    }
                    else
                    {
                        correctedEmissi = emissi*radiationCorrection;
                    }
                    
                }
                
            }
            
//             scalar T_Radiation = Foam::pow(GcInitial/(4.*Ste_Bol),0.25);
/*            
            if (T_Radiation>2000||T_Radiation<200)
            {
                Info<<" radiation abnormal: T "<<T_Radiation<<endl;
                Info<<" GcInitial: "<<GcInitial<<endl;
            }*/

            //h_coe = 0.0712114/(2.0*rb[3])*(2.0+1.1*0.890919*329.527*pow(2.0*rb[3],0.6));
            //h_coe = 112.737;from single particle code
//              h_coe = kappas/(this->d_)*(2.0*e_bed+0.69*pow(Res/e_bed,0.5)*cbrt(Prs)); //(2012 R. Mehrabian) its only for mass transfer

            label layerFlag = 0; //default case which is also the original case
            
            scalar layer3Fraction = mp_new[3] / this->mass0_;
            
            scalar layer2Fraction = mp_new[2] / this->mass0_;
            
            scalar dQTotal = 0;

            
            if (layer2Fraction < 0.01 && layer3Fraction  < 0.01)
            {
                layerFlag = 1; // ash and char layer dont resolve heat conduction
            }
            else if (layer3Fraction  < 0.01 && layer2Fraction > 0.01)
            {
                layerFlag = 2; // ash  layer dont resolve heat conduction
            }
            else if (layer3Fraction  > 0.01 && layer2Fraction < 0.005)
            {
                layerFlag = 4; // char layer dont resolve heat conduction, wery unlikely case
            }

            //Update all the boundary temperature at T+dt
// VID2_ = QComb_;
            //eq4 Tb[2], Qrxn,b2=0
            switch(layerFlag) 
            {
                case 1:
                    {
                        dQTotal = this->dT_ + QComb_;
                        
                        if (!usePpConduction)
                        {
                            dQTotal = QComb_;
                        }
//Info<<"layerFlag "<< layerFlag<<endl;           
// Info<<"this->dT_ "<< this->dT_<<endl;
// Info<<"QComb_ "<< QComb_<<endl;
// Info<<" TInifinit "<< TInifinit<<endl; 
// Info<<"Tb[3] "<<Tb[3] <<endl;
// Info<<"Tp[1] "<<Tp[1] <<endl;      
// Info<<"rb[3] "<<rb[3] <<endl;
// Info<<"rp[1] "<<rp[1] <<endl;  
                        if (particleShape == 1)
                        {
                            if (surfaceExplicitFlag)
                            {
                                bool haveRealRoots = eq4_cylinderL_check(h_coe, correctedEmissi, Ste_Bol, kp[1], Ab[3], rb[3], rp[1], TInifinit, Tp[1], GcInitial, dQTotal, Xi, Tb[3]);
                                
                                if(!haveRealRoots || lastIteration)
                                {
                                    Tb[3] = eq4_cylinderL_Explicit(h_coe, correctedEmissi, Ste_Bol, kp[1], Ab[3], rb[3], rp[1], TInifinit, Tp[1], Tb[3], GcInitial, dQTotal, Xi);
                                }
                            }
                            else
                            {
                                Tb[3] = eq4_cylinderL(h_coe, correctedEmissi, Ste_Bol, kp[1], Ab[3], rb[3], rp[1], TInifinit, Tp[1], GcInitial, dQTotal, Xi);
                            }
                        }
                        else
                        {
                            if (surfaceExplicitFlag)
                            {
                                bool haveRealRoots = eq4_check(h_coe, correctedEmissi, Ste_Bol, kp[1], Ab[3], rb[3], rp[1], TInifinit, Tp[1], GcInitial, dQTotal, Tb[3]);
                                
                                if(!haveRealRoots || lastIteration)
                                {
                                    Tb[3] = eq4_Explicit(h_coe, correctedEmissi, Ste_Bol, kp[1], Ab[3], rb[3], rp[1], TInifinit, Tp[1], Tb[3], GcInitial, dQTotal);
                                }
                            }
                            else
                            {
                                Tb[3] = eq4(h_coe, correctedEmissi, Ste_Bol, kp[1], Ab[3], rb[3], rp[1], TInifinit, Tp[1], GcInitial, dQTotal);
                            }
                        }
                    }
                    break;
                case 2:
                    {
                        dQTotal = this->dT_ + QComb_;
                        
                        if (!usePpConduction)
                        {
                            dQTotal = QComb_;
                        }
                        
                        if (particleShape == 1)
                        {
                            if (surfaceExplicitFlag)
                            {
                                bool haveRealRoots = eq4_cylinderL_check(h_coe, correctedEmissi, Ste_Bol, kp[2], Ab[3], rb[3], rp[2], TInifinit, Tp[2], GcInitial, dQTotal, Xi, Tb[3]);
                                
                                if(!haveRealRoots || lastIteration)
                                {
                                    Tb[3] = eq4_cylinderL_Explicit(h_coe, correctedEmissi, Ste_Bol, kp[2], Ab[3], rb[3], rp[2], TInifinit, Tp[2], Tb[3], GcInitial, dQTotal, Xi);
                                }
                            }
                            else
                            {
                                Tb[3] = eq4_cylinderL(h_coe, correctedEmissi, Ste_Bol, kp[2], Ab[3], rb[3], rp[2], TInifinit, Tp[2], GcInitial, dQTotal, Xi);
                            }
                        }
                        else
                        {
                            if (surfaceExplicitFlag)
                            {
                                bool haveRealRoots = eq4_check(h_coe, correctedEmissi, Ste_Bol, kp[2], Ab[3], rb[3], rp[2], TInifinit, Tp[2], GcInitial, dQTotal, Tb[3]);
                                
                                if(!haveRealRoots || lastIteration)
                                {
                                    Tb[3] = eq4_Explicit(h_coe, correctedEmissi, Ste_Bol, kp[2], Ab[3], rb[3], rp[2], TInifinit, Tp[2], Tb[3], GcInitial, dQTotal);
                                }
                            }
                            else
                            {
                                Tb[3] = eq4(h_coe, correctedEmissi, Ste_Bol, kp[2], Ab[3], rb[3], rp[2], TInifinit, Tp[2], GcInitial, dQTotal);
                            }
                        }
                        
                    }
                    break;
                case 4:
                    {
                        
                        dQTotal = this->dT_ + QComb_;
                        
                        if (!usePpConduction)
                        {
                            dQTotal = QComb_;
                        }
                        
                        if (particleShape == 1)
                        {
                            if (surfaceExplicitFlag)
                            {
                                bool haveRealRoots = eq4_cylinderL_check(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[2], TInifinit, Tp[2], GcInitial, dQTotal, Xi, Tb[3]);
                                
                                if(!haveRealRoots || lastIteration)
                                {
                                    Tb[3] = eq4_cylinderL_Explicit(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[2], TInifinit, Tp[2], Tb[3], GcInitial, dQTotal, Xi);
                                }
                            }
                            else
                            {
                                Tb[3] = eq4_cylinderL(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[2], TInifinit, Tp[2], GcInitial, dQTotal, Xi);
                            }
                        }
                        else
                        {
                            if (surfaceExplicitFlag)
                            {
                                bool haveRealRoots = eq4_check(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[2], TInifinit, Tp[2], GcInitial, dQTotal, Tb[3]);
                                
                                if(!haveRealRoots || lastIteration)
                                {
                                    Tb[3] = eq4_Explicit(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[2], TInifinit, Tp[2], Tb[3], GcInitial, dQTotal);
                                }
                            }
                            else
                            {
                                Tb[3] = eq4(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[2], TInifinit, Tp[2], GcInitial, dQTotal);
                            }
                        }
                    }
                    break;
                default:
                    {
                        dQTotal = this->dT_+ QComb_;
                        
                        if (!usePpConduction)
                        {
                            dQTotal = 0;
                        }
                        
                        if (particleShape == 1)
                        {
                            if (surfaceExplicitFlag)
                            {
                                bool haveRealRoots = eq4_cylinderL_check(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[3], TInifinit, Tp[3], GcInitial, dQTotal, Xi, Tb[3]);
                                
                                if(!haveRealRoots || lastIteration)
                                {
                                    Tb[3] = eq4_cylinderL_Explicit(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[3], TInifinit, Tp[3], Tb[3], GcInitial, dQTotal, Xi);
                                }
                            }
                            else
                            {
                                Tb[3] = eq4_cylinderL(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[3], TInifinit, Tp[3], GcInitial, dQTotal, Xi);
                            }
                        }
                        else
                        {
                            if (surfaceExplicitFlag)
                            {
                                bool haveRealRoots = eq4_check(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[3], TInifinit, Tp[3], GcInitial, dQTotal, Tb[3]);
                                
                                if(!haveRealRoots || lastIteration)
                                {
                                    Tb[3] = eq4_Explicit(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[3], TInifinit, Tp[3], Tb[3], GcInitial, dQTotal);
                                }
                            }
                            else
                            {
                                Tb[3] = eq4(h_coe, correctedEmissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[3], TInifinit, Tp[3], GcInitial, dQTotal);
                            }
                        }
                    }
            }
            
            //set a limit of maximum particle temperature changing avoid unreasonable tempearture jump !!!!!!!!!!!! need to be fixed TODO
            if (Tb[3] - Tb_old[3] > 50.0)
            {
                Tb[3] = Tb_old[3] + 50.0;
            }
            
            if (Tb[3]  > 1600.0 || Tb[3] < 200.0)
            {
                Tb[3] = 1473.0;
            }
            
            //debug info
//             if (Tb[3]>900)
//             {
//                 Info<<"iter_0 "<<iter_0<<endl; 
//                 Info<<"Tb[2] "<<Tb[2]<<endl;    
//                 Info<<"Tb[3] "<<Tb[3]<<endl; 
//                 Info<<"Tb_old[3] "<<Tb_old[3]<<endl; 
//                 Info<<"TInifinit "<<TInifinit<<endl; 
//                 Info<<"QComb "<<QComb_<<endl;
//                 Info<<"this->dT_ "<<this->dT_<<endl;
//                 Info<<"GcInitial "<<GcInitial<<endl;  
//                 Info<<"Tradia "<<pow(GcInitial/(4*Ste_Bol), 0.25)<<endl;  
//                 Info<<"contactedNeighbourNum "<<contactedNeighbourNum<<endl;  ppGc = 4*Ste_Bol*this->ppT4sum_/contactedNeighbourNum;
//                 Info<<"t_step "<<t_step<<endl; 
//                 Info<<"t_current "<<t_current<<endl;    
//                 Info<<"dt "<<dt<<endl; 
//             }
//             
            aa = correctedEmissi*Ste_Bol*Ab[3];
            VID1_ = h_coe*Ab[3]*(Tb[3] - TInifinit);
            VID2_ = (aa*Foam::pow4(Tb[3]) - aa*(GcInitial/(4.*Ste_Bol)));
            Qin_ = pow(GcInitial/(4.*Ste_Bol),0.25);
//             if(contactedNeighbourNum > 0.99)
//             {
//                 ppGc = (1.0 - e_bed)*4*Ste_Bol*this->ppT4sum_/contactedNeighbourNum + e_bed*GcInitial;
//                 
//             }
            
            VID4_ = this->dT_;
            
            QLayer_ = contactedNeighbourNum;
            switch(layerFlag) 
            {
                case 1:
                    {
                        VID3_  = d_dr_cylinderL(kp[1], Ab[3], rb[3], rp[1], Xi)*(Tb[3]-Tp[1]);
//                         Qin_  = d_dr_cylinderL(kp[1], Ab[3], rb[3], rp[1], Xi);
//                         QLayer_ = rp[1];
                    }
                break;
                case 2:
                    {
                        VID3_  = d_dr_cylinderL(kp[2], Ab[3], rb[3], rp[2], Xi)*(Tb[3]-Tp[2]);
//                         Qin_  = d_dr_cylinderL(kp[2], Ab[3], rb[3], rp[2], Xi);
//                         QLayer_ = rp[2];
                    }
                break;
                case 4:
                    {
                        VID3_  = d_dr_cylinderL(kp[3], Ab[3], rb[3], rp[2], Xi)*(Tb[3]-Tp[2]);
//                         Qin_  = d_dr_cylinderL(kp[3], Ab[3], rb[3], rp[2], Xi);
//                         QLayer_ = rp[2];
                    }
                break;
                default:
                    {
                        VID3_  = d_dr_cylinderL(kp[3], Ab[3], rb[3], rp[3], Xi)*(Tb[3]-Tp[3]);
//                         Qin_  = d_dr_cylinderL(kp[3], Ab[3], rb[3], rp[3], Xi);
//                         QLayer_ = rp[3];
                    }
            }
                    
            
            
//             VID1_ = this->ppT4sum_;
//             VID2_ = this->dT_;
//                 VID1_ = dHeatSRCarrier[4];
//                 VID2_ = dHeatSRCarrier[5];
//                 VID3_ = dHeatSRCarrier[6];
//                 VID4_ = dHeatSRCarrier[7];
//             if(lastIteration)
//             {
//                 Info<<"convection_VID1_: "<<VID1_<<endl;
//                 Info<<"conduction_VID3_: "<<VID3_<<endl;
//                 Info<<"radiation_VID2_: "<<VID2_<<endl;
//                 Info<<"Source_dQTotal: "<<dQTotal<<endl;
//             }
            
// if(Tb[3]<0)
// {
//     Info<<"Tb[3]: "<<Tb[3]<<endl;
// }
     //char conversion
            forAll(dMassSRSolid, i)
            {
                dMassSRSolid[i] = 0.0;
            }
            forAll(dMassSRCarrier, i)
            {
                dMassSRCarrier[i] = 0.0;
            }            
            forAll(dHeatSRCarrier, i)
            {
                dHeatSRCarrier[i] = 0.0;
            }   

            if (mp_old[2] <= this->mass0_*p_WoodLeft || noCombustion)
            {
                Tb[2] = Tp[3];
                rComb_ = 0.;
                QComb_ = 0.;
                        
            }
            else
            {
                
                switch(layerFlag) 
                {
                    case 1:
                        {
                            Tb[2] = Tb[3];
                        }
                        break;
                    case 2:
                        {
                            Tb[2] = Tb[3];
                        }
                        break;
                    case 4:
                        {
                            Tb[2] = Tp[3];
                        }
                        break;
                    default:
                        {
// Info<<"Tb2 equation:---------------"<<endl;
                            if (particleShape == 1)
                            {
                                Tb[2] = eq7_2_cylinderL(kp[2], kp[3], Ab[2], rb[2], rp[2], rp[3], Tp[2], Tp[3], 0, Xi);
                            }
                            else
                            {
                                Tb[2] = eq7_2(kp[2], kp[3], Ab[2], rb[2], rp[2], rp[3], Tp[2], Tp[3], 0);
                            }
                        }
                }

// if (Tb[2] > 1100.0 ||  Tb[2] < 300.0)
// {
//     Info<<"Tb[2] "<<Tb[2]<<endl;    
//     Info<<"Tb[3] "<<Tb[3]<<endl; 
//     Info<<"TInifinit "<<TInifinit<<endl; 
//     Info<<"QComb "<<QComb_<<endl;
//     Info<<"e_bed "<<e_bed<<endl;
// }
                //set a limit of maximum particle temperature changing avoid unreasonable tempearture jump !!!!!!!!!!!! need to be fixed TODO
                if (Tb[2] - Tb_old[2] > 100.0)
                {
                    Tb[2] = Tb_old[2] + 100;
                }
                
                if (Tb[2]  > 1600.0 || Tb[2] < 200.0)
                {
                    Tb[3] = 1473.0;
                }
//     if (Tb[3]>400)
//     {
//         Info<<"iter_0 "<<iter_0<<endl; 
//         Info<<"Tb[2] "<<Tb[2]<<endl;    
//         Info<<"Tb[3] "<<Tb[3]<<endl; 
//         Info<<"t_current "<<t_current<<endl;    
//     }
//             
                scalar currentDegree = (mp_new[0]*(1.0 - Ywater00) + mp_new[1])/(this->mass0_*(1.0 - Ywater00));
                
                if (limitedCombustionFlag)
                {
                    if (combustionCriterion < 0.0)
                    {
                        if (flagDevo_ == 1)
                        {
                            noCombustionDueToDevo = true;
                        }
                    }
                    else
                    {
                        if (currentDegree > combustionCriterion)
                        {
                            noCombustionDueToDevo = true;
                        }
                    }
                }
//                 VID2_ = Tb[2];
//                 VID3_ = Tb_old[2];
                // Calc mass and enthalpy transfer due to surface reactions  
                calcSurfaceReactions
                (
                    cloud,
                    td,
                    t_step,
                    rb[3],
                    rb[2],
                    Tb[2],
                    massit,
                    canCombust_,
                    Ne,
                    YMix,
                    YGas_,
                    YLiquid_,
                    YSolid_,
                    dMassSRGas,
                    dMassSRLiquid,
                    dMassSRSolid,
                    dMassSRCarrier,
                    Sh,
                    dhsTrans_dt,
                    QComb_,
                    Res,
                    TInifinit,
                    rhos,
                    mus,
                    Xi,
                    dHeatSRCarrier,
                    particleShape,
                    Csavg,
                    coarseGridSurfaceCombustion,
                    this->d_,
                    heatRatio,
                    noCombustionDueToDevo
                );
                
                rComb_ = mag(dMassSRSolid[idChar]);              
            }             
       
            //Devolatilisation
            forAll(dMassDV, i)
            {
                dMassDV[i] = 0.0;
            }
            forAll(dMassDVSoild, i)
            {
                dMassDVSoild[i] = 0.0;
            }            
            
            if (flagDevo_ == 1 )
            {
                switch(layerFlag) 
                {
                    case 1:
                        {
                             Tb[1] = Tb[2];
                        }
                        break;
                    case 2:
                        {
//                             Info<<"layerFlag 2 - Tb1 equation:---------------"<<endl;
                            if (particleShape == 1)
                            {
                                Tb[1] = eq7_2_cylinderL(kp[1], kp[2], Ab[1], rb[1], rp[1], rp[2], Tp[1], Tp[2], QDevo,Xi);
                            }
                            else
                            {
                                Tb[1] = eq7_2(kp[1], kp[2], Ab[1], rb[1], rp[1], rp[2], Tp[1], Tp[2], QDevo);
                            }
                        }
                        break;
                    case 4:
                        {
//                             Info<<"layerFlag 4 - Tb1 equation:---------------"<<endl;
                            if (particleShape == 1)
                            {
                                Tb[1] = eq7_2_cylinderL(kp[1], kp[2], Ab[1], rb[1], rp[1], rp[2], Tp[1], Tp[2], QDevo,Xi);
                            }
                            else
                            {
                                Tb[1] = eq7_2(kp[1], kp[2], Ab[1], rb[1], rp[1], rp[2], Tp[1], Tp[2], QDevo);
                            }
                        }
                        break;
                    default:
                        {
//                             Info<<"layerFlag default - Tb1 equation:---------------"<<endl;
                            if (particleShape == 1)
                            {
                                Tb[1] = eq7_2_cylinderL(kp[1], kp[2], Ab[1], rb[1], rp[1], rp[2], Tp[1], Tp[2], QDevo,Xi);
                            }
                            else
                            {
                                Tb[1] = eq7_2(kp[1], kp[2], Ab[1], rb[1], rp[1], rp[2], Tp[1], Tp[2], QDevo);
                            }
                
                        }
                }
                
                if (Tb[1] >= cloud.constProps().TDevol() && Tp[1] >= cloud.constProps().TDevol())
                {
                    deltRadiusSubLayer = (rb[1]-rb[0])/numberSubLayer;
                    
                    deltTemperature1 = (Tb[0]-Tp[1])/(rb[0]-rp[1]);
                    b0Temperature1 = (rb[0]*Tp[1]-rp[1]*Tb[0])/(rb[0]-rp[1]);
                    
                    deltTemperature2 = (Tp[1]-Tb[1])/(rp[1]-rb[1]);
                    b0Temperature2 = (rp[1]*Tb[1]-rb[1]*Tp[1])/(rp[1]-rb[1]);
                    
                    scalar massLayerCum = 0.0;
    /*Info<<"TInifinit "<<TInifinit <<endl;               
    Info<<"Tb[3] "<<Tb[3]<<endl;
    Info<<"Tb[2] "<<Tb[2]<<endl; 
    Info<<"Tb[1] "<<Tb[1]<<endl; */               
                    for (subLayer = 1; subLayer < int(numberSubLayer); subLayer++)
                    {
                        radiusSubLayerL = rb[0] + (subLayer-1)*deltRadiusSubLayer;
                        radiusSubLayerM = radiusSubLayerL + 0.5*deltRadiusSubLayer;
                        radiusSubLayerH = radiusSubLayerL + deltRadiusSubLayer;
                        forAll(dMassDVTemp, i)
                        {
                            dMassDVTemp[i] = 0.0;
                        }
                        forAll(dMassDVSoildTemp, i)
                        {
                            dMassDVSoildTemp[i] = 0.0;
                        }
                        
                        if (particleShape == 1)
                        {
                            massSubLayer = Vol_Rin_cylinderL(radiusSubLayerL, radiusSubLayerH,Xi)/Vp[1]*mp_new[1];
                        }
                        else
                        {
                            massSubLayer = Vol_Rin(radiusSubLayerL, radiusSubLayerH)/Vp[1]*mp_new[1];
                        }
                        
                        massLayerCum = massLayerCum + massSubLayer;
                        if (radiusSubLayerM<=rp[1])
                        {
                            temperatureSubLayer = radiusSubLayerM*deltTemperature1 + b0Temperature1;
                        }
                        else
                        {
                            temperatureSubLayer = radiusSubLayerM*deltTemperature2 + b0Temperature2;
                        }
                        
                        if (temperatureSubLayer > 1000.0)
                        {
                            temperatureSubLayer = 1000.0;
                        }

                        cloud.devolatilisation().calculate
                        (
                            t_step,
                            this->age_,
                            this->mass0_,
                            massSubLayer*YdryDB00,
                            temperatureSubLayer,
                            YMix[GAS]*YGas_,
                            YMix[LIQ]*YLiquid_,
                            YMix[SLD]*YSolid_, 
                            canCombust_,
                            dMassDVTemp,
                            dMassDVSoildTemp
                        );
                        dMassDV = dMassDV + dMassDVTemp;
                        dMassDVSoild = dMassDVSoild + dMassDVSoildTemp;
                    }
                    
                    forAll(dMassDVTemp, i)
                    {
                        dMassDVTemp[i] = 0.0;
                    }
                    forAll(dMassDVSoildTemp, i)
                    {
                        dMassDVSoildTemp[i] = 0.0;
                    }
            
                    massSubLayer = mp_new[1] - massLayerCum;
                    temperatureSubLayer =  Tb[1];
                    cloud.devolatilisation().calculate
                    (
                        t_step,
                        this->age_,
                        this->mass0_,
                        massSubLayer*YdryDB00,
                        temperatureSubLayer,
                        YMix[GAS]*YGas_,
                        YMix[LIQ]*YLiquid_,
                        YMix[SLD]*YSolid_, 
                        canCombust_,
                        dMassDVTemp,
                        dMassDVSoildTemp
                    );               
                    dMassDV = dMassDV + dMassDVTemp;
                    dMassDVSoild = dMassDVSoild + dMassDVSoildTemp;               
                    rDevo_ = mag(dMassDVSoild[idActiveDryWood]);
                    rChar_ = mag(dMassDVSoild[idChar]);
                    
                    QDevo = -rDevo_ * cloud.constProps().LDevol() / t_step;
//                     QLayer_ = QDevo;
                    
                    if (cloud.heatTransfer().BirdCorrection())
                    {

                        forAll(dMassDV, i)
                        {
                            const label id = composition.localToCarrierId(GAS, i);
                            const scalar W = composition.carrier().Wi(id);
                            const scalar Ni = dMassDV[i]/(Ab[3]*t_step*W);

                            Ne_temp += Ni;
                        }
                    }
                }
                else
                {
                    
                    dMassDV = 0.;
                    dMassDVSoild = 0.;
                    rDevo_ = 0.;
                    rChar_ = 0.;
                    QDevo  = 0.;
                }
                
                if (mp_old[1]+rDry_*ratioWoodMoist-rDevo_/YdryDB00<=this->mass0_*p_WoodLeft)
                {
                    flagDevo_ = 0;
                }                 
/*Info<<" rDevo_ "<< rDevo_<<endl;    
Info<<" rChar_ "<< rChar_<<endl; */           
            }
            else
            {
                Tb[1] = Tp[2];
                rDevo_ = 0.;
                rChar_ = 0.;
                QDevo  = 0.;
            }
// Info<<"rChar_ "<<rChar_<<endl;

            //Drying
            if (flagBoiling_ == 1)
            {
             //updat Fb, QDry_ (heat used to evaporate water )
                Fb1Value = Fb1(Tb[0]);
                if (Fb1Value>1.)
                {
                    Fb1Value = 1.;
                    Tb[0] = T_boiling;
                    if (particleShape == 1)
                    {
                        QDry_ = d_dr_cylinderL(kp[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1])-d_dr_cylinderL(kp[0],Ab[0],rb[0],rp[0],Xi)*(Tb[0]-Tp[0]);
                    }
                    else
                    {
                        QDry_ = d_dr(kp[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1])-d_dr(kp[0],Ab[0],rb[0],rp[0])*(Tb[0]-Tp[0]);
                    }
                    
                  
                    if (QDry_ < 0)
                    {
//                         FatalErrorInFunction
//                         <<"QDry_ < 0" << exit(FatalError);
                        QDry_ = 0;
                    }
                }
                else
                {
                    if (particleShape == 1)
                    {
                        Tb[0] = eq7_3_cylinderL(kp[0], kp[1], Ab[0], rb[0], rp[0], rp[1], Tp[0], Tp[1], Fb1Value,Xi);
                    }
                    else
                    {
                        Tb[0] = eq7_3(kp[0], kp[1], Ab[0], rb[0], rp[0], rp[1], Tp[0], Tp[1], Fb1Value);
                    }
                    
                    Fb1Value = Fb1(Tb[0]);
                                      
                    if (Tb[0] >= T_boiling)
                    {
                        Fb1Value = 1.;
                        Tb[0] = T_boiling;
                    }
                    
                    if (particleShape == 1)
                    {
                        QDry_ = d_dr_cylinderL(kp[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1])*Fb1Value;
                    }
                    else
                    {
                        QDry_ = d_dr(kp[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1])*Fb1Value;
                    }
                    
                    if (QDry_ < 0)
                    {
                        QDry_ = 0;
                    }
                }
        
                //Calculate drying rate
                rDry_ = QDry_*t_step/deltaHvap(Tb[0]);

                if (mp_old[0]-rDry_/Ywater00<=this->mass0_*p_WoodLeft)
                {
                    rDry_ = (mp_old[0]-this->mass0_*p_WoodLeft)*Ywater00;
                    flagBoiling_ = 0;
                } 
                
                if (cloud.heatTransfer().BirdCorrection())
                {
                    Ne_temp += rDry_/(Ab[3]*t_step*18.0);
                }
                
            }
            else
            {
                Tb[0] = Tp[1];
                QDry_ = 0.;
                Fb1Value = 0.;
                rDry_ = 0.;
            }
            
            Ne = Ne_temp;
            
            // mp0to1:mass transfer from wet wood to dry wood, mp1to2:mass transfer from dry wood to char
            mp0to1 = rDry_*ratioWoodMoist;                
            mp1to2 = rDevo_/YdryDB00*YashDB00+rChar_;
            mp2to3 = (ash_inchar_old/mp_old[2])*rComb_;
            // mGas0:gas released from wet wood, mGas1:gas released from dry wood
            mGas0 = rDry_;
            mGas1 = rDevo_/YdryDB00-mp1to2;
               
            // cpG00to1_temp, heat capacity of gas from wet wood, wet wood to dry wood
            // cpG01to2_temp, heat capacity of gas from wet wood, dry wood to char
            cpG00to1_temp = cp_water_vapor(0.5*(Tb[0]+Tb[1]));
            cpG01to2_temp = cp_water_vapor(0.5*(Tb[1]+Tb[2]));
            //   cpG02to3_temp = cp_water_vapor(0.5*(Tb[2]+Tb[3]));

            cpG11to2_temp = CP_GAS;
//             cpG12to3_temp = CP_GAS;
            
            // cpP0to1_temp, heat capacity of dry wood, wet wood to dry wood
            // cpP1to2_temp, heat capacity of char, dry wood to char
            cpP0to1_temp = cp_p_modified(1,0.5*(Tb[0]+Tp[1]),Ywater00);//T=(Tb[0]+Tp[1])/2
            cpP1to2_temp = cp_p_modified(2,0.5*(Tb[1]+Tb[2]),Ywater00);
            cpP2to3_temp = cp_p_modified(2,0.5*(Tb[2]+Tb[3]),Ywater00);
            
            // Qp0to1 Heat required when dry wood transfers from wet wood to dry wood
            // Qp1to2 Heat required when char transfers from dry wood to char
            Qp0to1 = mp0to1*cpP0to1_temp*(Tp[1]-Tb[0]);
            Qp1to2 = mp1to2*cpP1to2_temp*(Tp[2]-Tb[1]);
            Qp2to3 = mp2to3*cpP2to3_temp*(Tp[3]-Tb[2]);
            //              
            // Heat required when gas transfers from wet wood to drywood
            // Heat required when gas transfers from dry wood to char
            // Heat required when gas transfers from char to ash
            Qg0to1 = mGas0*cpG00to1_temp*(Tb[1]-Tb[0]);
            Qg1to2 = mGas0*cpG01to2_temp*(Tb[2]-Tb[1]) + mGas1*cpG11to2_temp*(Tb[2]-Tb[1]);
//             Qg2to3 = mGas0*cpG02to3_temp*(Tb[3]-Tb[2]) + mGas1*cpG12to3_temp*(Tb[3]-Tb[2]);            
            
            //Updat dTpdt at t+dt
            if (particleShape == 1)
            {
                switch(layerFlag) 
                {
                    case 1:
                        {
                            dTpdt[3] = 0;
                            //eq 9
                            dTpdt[2] = 0;
                            //eq 10
                            dTpdt[1] = (d_dr_cylinderL(alphap[1],Ab[1],rb[1],rp[1],Xi)*(Tb[1]-Tp[1])-d_dr_cylinderL(alphap[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1+QDevo*t_step)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        }
                        break;
                    case 2:
                        {
                            dTpdt[3] = 0;
                            //eq 9
                            dTpdt[2] = (d_dr_cylinderL(alphap[2],Ab[2],rb[2],rp[2],Xi)*(Tb[2]-Tp[2])-d_dr_cylinderL(alphap[2],Ab[1],rb[1],rp[2],Xi)*(Tb[1]-Tp[2]))/(Vp[2])-(Qg1to2+Qp1to2+QCombForLayer*t_step)/t_step/(rhop[2]*cp[2]*Vp[2]);
                            //eq 10
                            dTpdt[1] = (d_dr_cylinderL(alphap[1],Ab[1],rb[1],rp[1],Xi)*(Tb[1]-Tp[1])-d_dr_cylinderL(alphap[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        }
                        break;
                    case 4:
                        {
                            dTpdt[3] = (d_dr_cylinderL(alphap[3],Ab[3],rb[3],rp[3],Xi)*(Tb[3]-Tp[3])-d_dr_cylinderL(alphap[3],Ab[2],rb[2],rp[3],Xi)*(Tb[2]-Tp[3]))/(Vp[3])-(Qp2to3)/t_step/(rhop[3]*cp[3]*Vp[3]);
                            //eq 9
                            dTpdt[2] = (d_dr_cylinderL(alphap[2],Ab[2],rb[2],rp[2],Xi)*(Tb[2]-Tp[2])-d_dr_cylinderL(alphap[2],Ab[1],rb[1],rp[2],Xi)*(Tb[1]-Tp[2]))/(Vp[2])-(Qg1to2+Qp1to2)/t_step/(rhop[2]*cp[2]*Vp[2]);
                            //eq 10
                            dTpdt[1] = (d_dr_cylinderL(alphap[1],Ab[1],rb[1],rp[1],Xi)*(Tb[1]-Tp[1])-d_dr_cylinderL(alphap[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        }
                        break;
                    default:
                        {
                            dTpdt[3] = (d_dr_cylinderL(alphap[3],Ab[3],rb[3],rp[3],Xi)*(Tb[3]-Tp[3])-d_dr_cylinderL(alphap[3],Ab[2],rb[2],rp[3],Xi)*(Tb[2]-Tp[3]))/(Vp[3])-(Qp2to3)/t_step/(rhop[3]*cp[3]*Vp[3]);
                            //eq 9
                            dTpdt[2] = (d_dr_cylinderL(alphap[2],Ab[2],rb[2],rp[2],Xi)*(Tb[2]-Tp[2])-d_dr_cylinderL(alphap[2],Ab[1],rb[1],rp[2],Xi)*(Tb[1]-Tp[2]))/(Vp[2])-(Qg1to2+Qp1to2)/t_step/(rhop[2]*cp[2]*Vp[2]);
                            //eq 10
                            dTpdt[1] = (d_dr_cylinderL(alphap[1],Ab[1],rb[1],rp[1],Xi)*(Tb[1]-Tp[1])-d_dr_cylinderL(alphap[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]);
                             
                        }
                }
                
                //eq 11-1
                dTpdt[0] = d_dr_cylinderL((kp[0]/(cp[0]*rhop[0])),Ab[0],rb[0],rp[0],Xi)*(Tb[0]-Tp[0])/(Vp[0]);  
            }
            else
            {
                switch(layerFlag) 
                {
                    case 1:
                        {
                            dTpdt[3] = 0;
                            //eq 9
                            dTpdt[2] = 0;
                            //eq 10
                            dTpdt[1] = (d_dr(alphap[1],Ab[1],rb[1],rp[1])*(Tb[1]-Tp[1])-d_dr(alphap[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1+QDevo*t_step)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        }
                        break;
                    case 2:
                        {
                            dTpdt[3] = 0;
                            //eq 9
                            dTpdt[2] = (d_dr(alphap[2],Ab[2],rb[2],rp[2])*(Tb[2]-Tp[2])-d_dr(alphap[2],Ab[1],rb[1],rp[2])*(Tb[1]-Tp[2]))/(Vp[2])-(Qg1to2+Qp1to2+QCombForLayer*t_step)/t_step/(rhop[2]*cp[2]*Vp[2]);
                            //eq 10
                            dTpdt[1] = (d_dr(alphap[1],Ab[1],rb[1],rp[1])*(Tb[1]-Tp[1])-d_dr(alphap[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        }
                        break;
                    case 4:
                        {
                            dTpdt[3] = (d_dr(alphap[3],Ab[3],rb[3],rp[3])*(Tb[3]-Tp[3])-d_dr(alphap[3],Ab[2],rb[2],rp[3])*(Tb[2]-Tp[3]))/(Vp[3])-(Qp2to3)/t_step/(rhop[3]*cp[3]*Vp[3]);
                            //eq 9
                            dTpdt[2] = (d_dr(alphap[2],Ab[2],rb[2],rp[2])*(Tb[2]-Tp[2])-d_dr(alphap[2],Ab[1],rb[1],rp[2])*(Tb[1]-Tp[2]))/(Vp[2])-(Qg1to2+Qp1to2)/t_step/(rhop[2]*cp[2]*Vp[2]);
                            //eq 10
                            dTpdt[1] = (d_dr(alphap[1],Ab[1],rb[1],rp[1])*(Tb[1]-Tp[1])-d_dr(alphap[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        }
                        break;
                    default:
                        {
                            dTpdt[3] = (d_dr(alphap[3],Ab[3],rb[3],rp[3])*(Tb[3]-Tp[3])-d_dr(alphap[3],Ab[2],rb[2],rp[3])*(Tb[2]-Tp[3]))/(Vp[3])-(Qp2to3)/t_step/(rhop[3]*cp[3]*Vp[3]);
                            //eq 9
                            dTpdt[2] = (d_dr(alphap[2],Ab[2],rb[2],rp[2])*(Tb[2]-Tp[2])-d_dr(alphap[2],Ab[1],rb[1],rp[2])*(Tb[1]-Tp[2]))/(Vp[2])-(Qg1to2+Qp1to2)/t_step/(rhop[2]*cp[2]*Vp[2]);
                            //eq 10
                            dTpdt[1] = (d_dr(alphap[1],Ab[1],rb[1],rp[1])*(Tb[1]-Tp[1])-d_dr(alphap[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        }
                }
                //eq 11-1
                dTpdt[0] = d_dr((kp[0]/(cp[0]*rhop[0])),Ab[0],rb[0],rp[0])*(Tb[0]-Tp[0])/(Vp[0]); 
                
            }
                     
            switch(layerFlag) 
            {
                case 1:
                    {
                        for ( i=2;i<=3; i++)
                        {
                            Tp_lastiter[i]=Tp[i];
                            Tp[i] = Tb[i];       
                            tolerance_Tp[i] = VSMALL;
                            tolerance_Tb[i] = mag((Tb_lastiter[i] - Tb[i])/Tb[i]); 
                        } 
                    }
                    break;
                case 2:
                    {
                        Tp_lastiter[3]=Tp[3];
                        Tp[3] = Tb[3];       
                        tolerance_Tp[3] = VSMALL;
                        tolerance_Tb[3] = mag((Tb_lastiter[3] - Tb[3])/Tb[3]); 
                        Tp_lastiter[2]=Tp[2];
                        Tp[2] = Tp_old[2] + dTpdt[2]*t_step;;       
                        tolerance_Tp[2] = mag((Tp_lastiter[2] - Tp[2])/Tp[2]);
                        tolerance_Tb[2] = mag((Tb_lastiter[2] - Tb[2])/Tb[2]); 
                    }
                    break;
                case 4:
                    {
                        for ( i=2;i<=3; i++)
                        {
                            Tp_lastiter[i]=Tp[i];
                            Tp[i] = Tp_old[i] + dTpdt[i]*t_step;       
                            tolerance_Tp[i] = mag((Tp_lastiter[i] - Tp[i])/Tp[i]);
                            tolerance_Tb[i] = mag((Tb_lastiter[i] - Tb[i])/Tb[i]); 
                        }  
                    }
                    break;
                default:
                    {
                        for ( i=2;i<=3; i++)
                        {
                            Tp_lastiter[i]=Tp[i];
                            Tp[i] = Tp_old[i] + dTpdt[i]*t_step;       
                            tolerance_Tp[i] = mag((Tp_lastiter[i] - Tp[i])/Tp[i]);
                            tolerance_Tb[i] = mag((Tb_lastiter[i] - Tb[i])/Tb[i]); 
                        }  
                    }
                }
                
                if (Tp[2]<273.0 && layerFlag == 2)
                {
                    Tp[2] = Tp_old[2];
                    QCombForLayer = QComb_;
                    QComb_ = 0.0;
                }
                else
                {
                    QCombForLayer = 0.0;
                }
// if (Tb[3]  < 297)
// {
// //     Info<<"h_coe "<<h_coe<<endl;
//     Info<<"layerFlag "<< layerFlag<<endl;       
//     Info<<"Tb[3] "<<Tb[3] <<endl;
//     Info<<"Tp[3] "<<Tp[3] <<endl;
//     Info<<"Tb[2] "<<Tb[2] <<endl;
//     Info<<"Tp[2] "<<Tp[2] <<endl;
//     Info<<"Tb[1] "<<Tb[1] <<endl;
//     Info<<"Tp[1] "<<Tp[1] <<endl;
//     Info<<"Tb[0] "<<Tb[0] <<endl;
//     Info<<"Tp[0] "<<Tp[0] <<endl;
//     
//     Info<<"rhop[3]  "<<rhop[3]  <<endl;
//     Info<<"rhop[2]  "<<rhop[2]  <<endl;
// 
// //     Info<<"QComb_ "<<QComb_<<endl;
// //     Info<<"rComb_ "<<rComb_<<endl;    
// //     Info<<"noCombustionDueToDevo "<<noCombustionDueToDevo<<endl;
// //     
// //     Info<<" TInifinit "<< TInifinit<<endl; 
// //     Info<<"Gc "<<Gc <<endl;
// //     Info<<"dQTotal "<<dQTotal<<endl;
// //     
// //     Info<<"flagBoiling_ "<<flagBoiling_<<endl;
// //     Info<<"flagDevo_ "<<flagDevo_ <<endl;
// //     Info<<"QDevo "<<QDevo<<endl;
// //     Info<<"rDevo_  "<<rDevo_ <<endl;
//     Info<<"dTpdt[3]  "<<dTpdt[3]  <<endl;
//     Info<<"dTpdt[2]  "<<dTpdt[2]  <<endl;
//     /*Info<<"rb[2]  "<<rb[2]  <<endl;
//     Info<<"rp[2]  "<<rp[2]  <<endl;
//     Info<<"rb[1]  "<<rb[1]  <<endl;
//     Info<<"alphap[2]  "<<alphap[2]  <<endl;
//     Info<<"kp[2]  "<<kp[2]  <<endl;
//     Info<<"rhop[2]  "<<rhop[2]  <<endl;
//     Info<<"cp[2]  "<<cp[2]  <<endl;
//     Info<<"Ab[2]  "<<Ab[2]  <<endl;
//     Info<<"Tb[2]-Tp[2]  "<<Tb[2]-Tp[2]  <<endl;
//     Info<<"d_dr(alphap[2],Ab[2],rb[2],rp[2])  "<<d_dr(alphap[2],Ab[2],rb[2],rp[2])  <<endl;
//     Info<<"(d_dr(alphap[2],Ab[2],rb[2],rp[2])*(Tb[2]-Tp[2]))/(Vp[2])  "<<(d_dr(alphap[2],Ab[2],rb[2],rp[2])*(Tb[2]-Tp[2]))/(Vp[2])  <<endl;
//     Info<<"Tb[1]-Tp[2]  "<<Tb[1]-Tp[2]  <<endl;
//     Info<<"-d_dr(alphap[2],Ab[1],rb[1],rp[2])*(Tb[1]-Tp[2])/(Vp[2])  "<<-d_dr(alphap[2],Ab[1],rb[1],rp[2])*(Tb[1]-Tp[2])/(Vp[2])  <<endl;
//     Info<<"(Qg1to2+Qp1to2)/t_step/(rhop[2]*cp[2]*Vp[2]) "<<(Qg1to2+Qp1to2)/t_step/(rhop[2]*cp[2]*Vp[2]) <<endl;
//     Info<<"dTpdt[1]  "<<dTpdt[1] <<endl;
//     
//     Info<<"(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]) "<<(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]) <<endl;
//     Info<<"dTpdt[0] "<<dTpdt[0]  <<endl;
//    */ 
// 
//     
//     Info<<endl;
// }


            //update Tp and tolerance of Tp and Tb 
            for ( i=0;i<=1; i++)
            {
                Tp_lastiter[i]=Tp[i];
                Tp[i] = Tp_old[i] + dTpdt[i]*t_step;       
                tolerance_Tp[i] = mag((Tp_lastiter[i] - Tp[i])/Tp[i]);
                tolerance_Tb[i] = mag((Tb_lastiter[i] - Tb[i])/Tb[i]);  
                
            }   
            if (flagBoiling_ == 0)
            {
                Tb[0]=Tp[1];
                Tp[0]=Tp[1];
                tolerance_Tp[0] = VSMALL;
                tolerance_Tb[0] = VSMALL;
            }
            if (flagDevo_ == 0)
            {
                Tb[1]=Tp[2];
                Tp[1]=Tp[2];
                tolerance_Tp[1] = VSMALL;
                tolerance_Tb[1] = VSMALL;
            }
            if (mp_old[2] <= this->mass0_*p_WoodLeft)
            {
                Tb[2]=Tp[3];
                Tp[2]=Tp[3];
                tolerance_Tb[2] = VSMALL;
                tolerance_Tp[2] = VSMALL;
            }
            
            if (flagBoiling_ == 1  && (mp_old[1]+rDry_*ratioWoodMoist-rDevo_/YdryDB00 > this->mass0_*p_WoodLeft))
            {
                flagDevo_ = 1;
            }   
            
            
// if (iter_0>50)
// {
//     for ( i=0;i<=3; i++)
//             {
//                 
//                 Info<<"tolerance_Tp["<<i<<"]: "<< mag((Tp_lastiter[i] - Tp[i])/Tp[i])<<"   deltTolerance: "<<tolerance_Tp[i] - tolerance<<"  Tp["<<i<<"]: "<<Tp[i]<<"  Tp_lastiter["<<i<<"]: "<<Tp_lastiter[i]<<endl;
//                 Info<<"tolerance_Tb["<<i<<"]: "<< mag((Tb_lastiter[i] - Tb[i])/Tb[i])<<"   deltTolerance: "<<tolerance_Tb[i] - tolerance<<"  Tb["<<i<<"]: "<<Tb[i]<<"  Tb_lastiter["<<i<<"]: "<<Tb_lastiter[i]<<endl;
//                 
//             }    
// }
            
            iter_0 = iter_0+1;
            
            if ((iter_0+1) == maxIters)
            {
                lastIteration = true;
            }
        }while (((tolerance_Tb[0] > tolerance || tolerance_Tb[1] > tolerance || tolerance_Tb[2] > tolerance
        || tolerance_Tb[3] > tolerance  || tolerance_Tp[0] > tolerance || tolerance_Tp[1] > tolerance 
        || tolerance_Tp[2] > tolerance  || tolerance_Tp[3] > tolerance)  && iter_0 < maxIters) || iter_0<2);

// if (Tb[2]>550)
// {
//     Info<<"Tb[2]: "<<Tb[2]<<endl;
//     Info<<"heat of C + O2: "<<dHeatSRCarrier[0]<<endl;
//     Info<<"heat of C + CO2: "<<dHeatSRCarrier[1]<<endl;
//     Info<<"heat of C + H2O: "<<dHeatSRCarrier[2]<<endl;
//     Info<<"heat of C + H2: "<<dHeatSRCarrier[3]<<endl;
// }
        
        //- Layer mass, radius, volume
        mp[0] = mp_old[0]-rDry_/Ywater00; 
        mp[1] = mp_old[1]+rDry_*ratioWoodMoist-rDevo_/YdryDB00;
        mp[2] = mp_old[2]+rDevo_/YdryDB00*YashDB00+rChar_-rComb_;
        ash_inchar = ash_inchar_old + rDevo_/YdryDB00*YashDB00 - (ash_inchar_old/mp_old[2])*rComb_;
        mp[3] = mp_old[3]+(ash_inchar_old/mp_old[2])*rComb_; 
        
        Vp[0] = (Vp_old[0]/mp_old[0])*mp[0];
        Vp1_temp = (Vp_old[1]/mp_old[1])*(mp_old[1]-rDevo_/YdryDB00);
        Vp[1] = Vp1_temp+(1. - drySrin)*(Vp_old[0] - Vp[0]);
        Vp2_temp = (Vp_old[2]/mp_old[2])*(mp_old[2]-rComb_);
        Vp[2] = Vp2_temp + (1. - devoSrin)*(Vp_old[1] - Vp1_temp);
        Vp[3] = Vp_old[3] + (1. - charSrin)*(Vp_old[2] - Vp2_temp);
        
        if (particleShape == 1)
        {
            rb[0] = radiusForCylinder(Xi,Vp[0]);
            rp[0] = R_Par_cylinderL(0.0,rb[0],Xi);
            Ab[0] = Area_cylinderL(rb[0],Xi);
            rb[1] = radiusForCylinder(Xi,Vp[0]+Vp[1]);
            rp[1] = R_Par_cylinderL(rb[0],rb[1],Xi);
            Ab[1] = Area_cylinderL(rb[1],Xi);
            rb[2] = radiusForCylinder(Xi,Vp[0]+Vp[1]+Vp[2]);
            rp[2] = R_Par_cylinderL(rb[1],rb[2],Xi);
            Ab[2] = Area_cylinderL(rb[2],Xi);
            rb[3] = radiusForCylinder(Xi,Vp[0]+Vp[1]+Vp[2]+Vp[3]);
            rp[3] = R_Par_cylinderL(rb[2],rb[3],Xi);
            Ab[3] = Area_cylinderL(rb[3],Xi);
        }
        else
        {
            rb[0] = Foam::cbrt(0.75*(Vp[0])/constant::mathematical::pi);
            rp[0] = R_Par(0.0,rb[0]);
            Ab[0] = Area_Sph(rb[0]);
            rb[1] = Foam::cbrt(0.75*(Vp[0]+Vp[1])/constant::mathematical::pi);
            rp[1] = R_Par(rb[0],rb[1]);
            Ab[1] = Area_Sph(rb[1]);
            rb[2] = Foam::cbrt(0.75*(Vp[0]+Vp[1]+Vp[2])/constant::mathematical::pi);
            rp[2] = R_Par(rb[1],rb[2]);
            Ab[2] = Area_Sph(rb[2]);
            rb[3] = Foam::cbrt(0.75*(Vp[0]+Vp[1]+Vp[2]+Vp[3])/constant::mathematical::pi);
            rp[3] = R_Par(rb[2],rb[3]);
            Ab[3] = Area_Sph(rb[3]); 
        
        }

        //- Sh update  Sh is not needed in this code
        //Sh = Sh - QDry_/t_step;//sign
        //Sh = Sh - rDevo_*cloud.constProps().LDevol()/t_step;
        //Sh = Sh - QComb_/t_step;// may not be correct, may better calculate it in surfacreaction function

        //- Mass and fraction update
        
        //drying
        dMassPC[idwater] = rDry_;
        dMassPCSoild[idDryWood] = rDry_*ratioWoodMoist;
        dMassPCSoild[idActiveDryWood] = -rDry_*ratioWoodMoist;         
        dMassLiquid = dMassLiquid + dMassPC;      
        massit = updateMassFractions(massit, dMassGasZero, dMassPC, dMassPCSoild);
        
        //Devolatilisation
        dMassGas = dMassGas + dMassDV;
        massit = updateMassFractions(massit, dMassGasZero, dMassLiquidZero, dMassDVSoild);

        //Combustion
        dMassGas = dMassGas + dMassSRGas;
        dMassLiquid = dMassLiquid + dMassSRLiquid;
        dMassSolid = dMassSolid + dMassSRSolid;
        massit = updateMassFractions(massit, dMassGasZero, dMassLiquidZero, dMassSolid);
        dMassSRCarrierTot = dMassSRCarrierTot + dMassSRCarrier;
        
        //- Add to cumulative phase change mass
        cloud.phaseChange().addToPhaseChangeMass(np0*rDry_);//need to be corrected
        cloud.devolatilisation().addToDevolatilisationMass(np0*rDevo_);

        //- Update molar emissions//Ts????
//         if (cloud.heatTransfer().BirdCorrection())
//         {
//             // Average molecular weight of carrier mix - assumes perfect gas
//             const scalar Wc = this->rhoc_*RR*this->Tc_/pc;
//             
//             forAll(dMassPC, i)
//             {
//                 const label idc = composition.localToCarrierId(idL, i);
//                 //                 const label idl = composition.globalIds(idL)[i];//idl to i?? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!??
//                 
//                 const scalar Cp = composition.carrier().Cp(idc, pc, Ts);
//                 const scalar W = composition.carrier().Wi(idc);
//                 const scalar Ni = dMassPC[i]/(this->areaS(2*rb[2])*t_step*W);
//                 
//                 const scalar Dab = composition.liquids().properties()[i].D(pc, Ts, Wc);
//                 
//                 // Molar flux of species coming from the particle (kmol/m^2/s)
//                 Ne += Ni;
//                 
//                 // Sum of Ni*Cpi*Wi of emission species
//                 NCpW += Ni*Cp*W;
//                 
//                 // Concentrations of emission species
//                 Cs[idc] += Ni*(2*rb[2])/(2.0*Dab);
//             }
//             
//             
//             // Note: hardcoded gaseous diffusivities for now
//             // TODO: add to carrier thermo
//             const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));
//             
//             forAll(dMassDV, i)
//             {
//                 const label id = composition.localToCarrierId(idG, i);
//                 const scalar Cp = composition.carrier().Cp(id, pc, Ts);
//                 const scalar W = composition.carrier().Wi(id);
//                 const scalar Ni = dMassDV[i]/(this->areaS(2*rb[2])*t_step*W);
//                 
//                 // Dab calc'd using API vapour mass diffusivity function
//                 const scalar Dab =
//                 3.6059e-3*(pow(1.8*Ts, 1.75))
//                 *sqrt(1.0/W + 1.0/Wc)
//                 /(pc*beta);
//                 
//                 Ne += Ni;
//                 NCpW += Ni*Cp*W;
//                 Cs[id] += Ni*(2*rb[2])/(2.0*Dab);
//             }      
//         }        
        //- Heat transfer
        //- update Sph

        Sph += t_step*h_coe*Ab[3];

        if (particleToGasOneWayHtc )
        {
            dhsTrans_dt = 0.0;
        }
        else
        {
//             dhsTrans_dt +=  t_step*(h_coe*Ab[3]*((Tb[3]+Tb_old[3])/2.0 - TInifinit) + (aa*Foam::pow4((Tb[3]+Tb_old[3])/2.0) - aa*(e_bed*GcInitial/(4.*Ste_Bol))));
            dhsTrans_dt +=  t_step*h_coe*Ab[3]*((Tb[3]+Tb_old[3])/2.0 - TInifinit);
        }
        
//         dhsTrans_dt +=  t_step*h_coe*Ab[3]*((Tb[3]+Tb_old[3])/2.0 - TInifinit);
        
        dhsTrans = dhsTrans + dhsTrans_dt;
        
//         QLayer_ += dhsTrans_dt;
// VID1_ = (Tb[3]+Tb_old[3])/2.0;
// VID2_ = TInifinit;
//         VID1_ = h_coe*Ab[3]*((Tb[3]+Tb_old[3])/2.0 - TInifinit);
//         VID4_ = dhsTrans;
        
        //- Update time
//debug info
      

// VID1_ += dhsTrans_dt;
        t_current = t_current + t_step;
        cumTime_ = cumTime_ + t_step;
            
        if (t_current + t_step >= dt)
        {
            t_step = dt - t_current;
        }
// Info<<"t_current out itr loop: "<<t_current<<endl;  
// Info<<"dt out itr loop: "<<dt<<endl;  
        this->T_ = /*Tb[2]*/Tb[3];
    }while (t_current<dt);
// Info<<"particle loop finished"<<endl;  
// Info<<"td.pc: "<<td.pc()<<endl;
// Info<<"td.Tc: "<<td.Tc()<<endl;
// Info<<"td.Uc: "<<td.Uc()<<endl;

    
    // equivalent diameter for cylinder and sphericity
    if (particleShape == 1)
    {
//         const scalar aspectRatio = (2*rb[3] - Xi)/(2*rb[3]);
        equivalent_d = Foam::cbrt(6*(Vp[0]+Vp[1]+Vp[2]+Vp[3])/constant::mathematical::pi);
//         scalar sphericity = (constant::mathematical::pi*equivalent_d*equivalent_d)/Ab[3];
        this->d_ = equivalent_d;
    }
    else
    {
        this->d_ = 2*rb[3]/*rb[2]*/;
    }
    //- Update class value
    this->Tp0_ = Tp[0];
    this->Tp1_ = Tp[1];
    this->Tp2_ = Tp[2];
    this->Tp3_ = Tp[3];
    this->Tb0_ = Tb[0];
    this->Tb1_ = Tb[1];
    this->Tb2_ = Tb[2];
    this->Tb3_ = Tb[3];
    this->rb0_ = rb[0];
    this->rb1_ = rb[1];
    this->rb2_ = rb[2];
    this->rb3_ = rb[3];
    this->mp0_ = mp[0];
    this->mp1_ = mp[1];
    this->mp2_ = mp[2];
    this->mp3_ = mp[3];
    this->T_ = /*Tb[2]*/Tb[3];
    
    this->rho_ = (mp[0]+mp[1]+mp[2]+mp[3])/(Vp[0]+Vp[1]+Vp[2]+Vp[3]);
    this->ash_inchar_t_ = ash_inchar;

    //- Correct surface values due to emitted species
    //no correctSurfaceValues if BirdCorrection is not actived, for averaged method, this step is not considered.
    // cellValueSourceCorrection will be called in the move function in kinetic parcel part. here for averaged method this step is not implemented.
    if (coarseGrid)
    {
            //correct Reynolds number
        if (UScheme == "cellPoint")
        {
            Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
        }
        else
        {
            Res = this->Re(rhos, U0, cloud.Uavg()[this->cell()], this->d_, mus);
        }
            
    }
    else
    {
        // correct surface values
        this->correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);   
        
    }
    
    //- Calculate new particle velocity
    

   
//     // Remove the particle when mass falls below minimum threshold
//     if (np0*mass1 < cloud.constProps().minParcelMass())
//     {
//         td.keepParticle = false;
// 
//         if (cloud.solution().coupled())
//         {
//             scalar dm = np0*mass0;
// 
//             // Absorb parcel into carrier phase
//             forAll(YGas_, i)
//             {
//                 label gid = composition.localToCarrierId(GAS, i);
//                 cloud.rhoTrans(gid)[this->cell()] += dm*YMix[GAS]*YGas_[i];
//             }
//             forAll(YLiquid_, i)
//             {
//                 label gid = composition.localToCarrierId(LIQ, i);
//                 cloud.rhoTrans(gid)[this->cell()] += dm*YMix[LIQ]*YLiquid_[i];
//             }
// 
//             // No mapping between solid components and carrier phase
//             /*
//             forAll(YSolid_, i)
//             {
//                 label gid = composition.localToCarrierId(SLD, i);
//                 cloud.rhoTrans(gid)[this->cell()] += dm*YMix[SLD]*YSolid_[i];
//             }
//             */
// 
//             cloud.UTrans()[this->cell()] += dm*U0;
// 
//             cloud.hsTrans()[this->cell()] +=
//                 dm*HsEff(cloud, td, pc, T0, idG, idL, idS);
// 
//             cloud.phaseChange().addToPhaseChangeMass(np0*mass1);
//         }
// 
//         return;
//     }


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (cloud.solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(YGas_, i)
        {
            scalar dm = np0*dMassGas[i]; //in fact this is from devo in this code
            label gid = composition.localToCarrierId(GAS, i);
            scalar hs = composition.carrier().Hs(gid, pc, Tb2_); //Tp may cause some error
            if (particleToGasOneWayHtc)
            {
                dm = 0.0;
                hs = 0.0;
            }
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
//             VID4_  += dm*hs;
//             Info<<"Ygas: "<<endl;
//             Info<<"gid: "<<gid<<";  dm: "<<dm<<";  hs: "<<hs<<endl;
        }
        forAll(dMassSRCarrierTot, i) //Combustion
        {
            scalar dm = np0*dMassSRCarrierTot[i];
            scalar hs = composition.carrier().Hs(i, pc, Tb2_); //Tp may cause some error
            cloud.rhoTrans(i)[this->cell()] += dm;
//             if (particleToGasOneWayHtc)
//             {
//                 dm = 0.0;
//                 hs = 0.0;
//             }
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
//             VID4_  += dm*hs;
//             Info<<"dMassSRCarrierTot: "<<endl;
//             Info<<"i: "<<i<<";  dm: "<<dm<<";  hs: "<<hs<<endl;
        }
        forAll(YLiquid_, i)
        {
            scalar dm = np0*dMassLiquid[i]; //in fact this is from drying in this code
            label gid = composition.localToCarrierId(LIQ, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            if (particleToGasOneWayHtc)
            {
                dm = 0.0;
                hs = 0.0;
            }
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
//             VID4_  += dm*hs;
//             Info<<"YLiquid_: "<<endl;
//             Info<<"gid: "<<gid<<";  dm: "<<dm<<";  hs: "<<hs<<endl;
//             scalar ha = composition.carrier().Ha(gid, pc, T0);
//             Info<<"gid: "<<gid<<";  ha: "<<ha<<endl;
        }

        // Update momentum transfero
        cloud.UTrans()[this->cell()] += np0*dUTrans;
//         cloud.UTrans()[this->cell()] += Zero;
        cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTrans()[this->cell()] += np0*dhsTrans;
        cloud.hsCoeff()[this->cell()] += np0*Sph;
//         VID4_  += np0*dhsTrans;
        // Update radiation fields
        if (cloud.radiation())
        {
            scalar ap;
            if (particleShape == 1)
            {
                ap = constant::mathematical::pi*rb[3]*(3*rb[3]-Xi)/2;;
            }
            else
            {
                ap = this->areaP();
            }
            const scalar T4 = pow4(Tb3_);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;
        }
    }
}

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcV
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const bool coarseGrid = cloud.constProps().coarseGrid();
    
    word UScheme = cloud.solution().interpolationSchemes().lookup("U");
    word TScheme = cloud.solution().interpolationSchemes().lookup("T");
    
    scalarField mp(4);
    
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    mp[0] = this->mp0_;
    mp[1] = this->mp1_;
    mp[2] = this->mp2_;
    mp[3] = this->mp3_;
    
    scalar TInifinit = td.Tc();
    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;
    
    scalar Ts, rhos, mus, Prs, kappas;
    
    scalar Res;
    
     if (coarseGrid)
        {
            // Calc surface values
            if (TScheme == "cellPoint")
            {
                TInifinit = cloud.T()[this->cell()];
            }
            else
            {
                TInifinit = cloud.Tavg()[this->cell()];
            }

            // Surface temperature using two thirds rule
            Ts = (2.0*T0 + TInifinit)/3.0;

            if (Ts < cloud.constProps().TMin())
            {
                if (debug)
                {
                    WarningInFunction
                        << "Limiting parcel surface temperature to "
                        << cloud.constProps().TMin() <<  nl << endl;
                }

                Ts = cloud.constProps().TMin();
            }

            // Assuming thermo props vary linearly with T for small d(T)
            const scalar TRatio = TInifinit/Ts;

            rhos = cloud.rhoavg()[this->cell()]*TRatio;

            mus = cloud.muavg()[this->cell()];;
            
            // Reynolds number
            if (UScheme == "cellPoint")
            {
                Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
            }
            else
            {
                Res = this->Re(rhos, U0, cloud.Uavg()[this->cell()], this->d_, mus);
            }
            
            
        }
        else
        {
            // Calc surface values
            this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
            
            // Reynolds number
            Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
            
        }
    
    
    this->U_ = this->calcVelocity(cloud, td, dt, Res, mus, mp[0]+mp[1]+mp[2]+mp[3], Su, dUTrans, Spu);
}
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcDevolatilisation
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar age,
    const scalar Ts,
    const scalar d,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV,
    scalarField& dMassDVSoild,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
) const
{
    // Check that model is active
    if (!cloud.devolatilisation().active())
    {
        if (canCombust != -1)
        {
            canCombust = 1;
        }
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().TDevol();
    (void)cloud.constProps().LDevol();

    // Check that the parcel temperature is within necessary limits for
    // devolatilisation to occur
    if (T < cloud.constProps().TDevol() || canCombust == -1)
    {
        return;
    }

    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();


    // Total mass of volatiles evolved
    cloud.devolatilisation().calculate
    (
        dt,
        age,
        mass0,
        mass,
        T,
        YGasEff,
        YLiquidEff,
        YSolidEff,
        canCombust,
        dMassDV,
        dMassDVSoild
    );

    scalar dMassTot = sum(dMassDV);

    cloud.devolatilisation().addToDevolatilisationMass
    (
        this->nParticle_*dMassTot
    );

    Sh -= dMassTot*cloud.constProps().LDevol()/dt;

    // Update molar emissions
    if (cloud.heatTransfer().BirdCorrection())
    {
        // Molar average molecular weight of carrier mix
        const scalar Wc = max(small, td.rhoc()*RR*td.Tc()/td.pc());

        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));

        forAll(dMassDV, i)
        {
            const label id = composition.localToCarrierId(GAS, i);
            const scalar Cp = composition.carrier().Cp(id, td.pc(), Ts);
            const scalar W = composition.carrier().Wi(id);
            const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);

            // Dab calc'd using API vapour mass diffusivity function
            const scalar Dab =
                3.6059e-3*(pow(1.8*Ts, 1.75))
               *sqrt(1.0/W + 1.0/Wc)
               /(td.pc()*beta);

            N += Ni;
            NCpW += Ni*Cp*W;
            Cs[id] += Ni*d/(2.0*Dab);
        }
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar d,
    const scalar di,
    const scalar T,
    const scalar mass,
    const label canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans,
    scalar& QComb,
    const scalar Re,
    const scalar Tc,
    const scalar rhoc,
    const scalar muc,
    const scalar Xi,
    scalarField& dHeatSRCarrier,
    const label particleShape,
    const scalarField& Csavg,
    const bool coarseGrid,
    const scalar deq,
    const bool heatRatio,
    const bool limitedCombustion
) const
{
//     Info<<" surfacereaction called here"<<endl;
    // Check that model is active
    if (!cloud.surfaceReaction().active())
    {
        return;
    }
    
    // no combustion if devo is still strong
    if (limitedCombustion)
    {
        dhsTrans = 0.0;
        QComb = 0.0;
        return;
    }
    
    
    // for steability set a minimum char reaction temperature at 600K
    if (T < 550.0)
    {
        dhsTrans = 0.0;
        QComb = 0.0;
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().hRetentionCoeff();
    (void)cloud.constProps().TMax();

    // Check that model is active
    if (canCombust != 1)
    {
        return;
    }
    
    scalar hReaction;
    
    if (coarseGrid)
    {
        // Update surface reactions
        hReaction = cloud.surfaceReaction().calculate
        (
            dt,
            this->cell(),
            d,
            T,
            Tc,
            td.pc(),
            rhoc,
            mass,
            YGas,
            YLiquid,
            YSolid,
            YMix,
            N,
            dMassSRGas,
            dMassSRLiquid,
            dMassSRSolid,
            dMassSRCarrier,
            di,
            muc,
            Re,
            Xi,
            dHeatSRCarrier,
            particleShape,
            Csavg,
            coarseGrid,
            deq
        );
    }
    else
    {
        // Update surface reactions
        hReaction = cloud.surfaceReaction().calculate
        (
            dt,
            this->cell(),
            d,
            T,
            td.Tc(),
            td.pc(),
            td.rhoc(),
            mass,
            YGas,
            YLiquid,
            YSolid,
            YMix,
            N,
            dMassSRGas,
            dMassSRLiquid,
            dMassSRSolid,
            dMassSRCarrier,
            di,
            td.muc(),
            Re,
            Xi,
            dHeatSRCarrier,
            particleShape,
            Csavg,
            coarseGrid,
            deq
        );
    }

    cloud.surfaceReaction().addToSurfaceReactionMass
    (
        this->nParticle_
       *(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    );

    if (heatRatio && T>1500.0)
    {
        dhsTrans = dHeatSRCarrier[0] - hReaction;
        QComb = -(dhsTrans+dHeatSRCarrier[1]+dHeatSRCarrier[2]+dHeatSRCarrier[3])/dt;
        
    }
    else if (heatRatio && T>1550.0)
    {
        dhsTrans = dHeatSRCarrier[0]+dHeatSRCarrier[1]+dHeatSRCarrier[2]+dHeatSRCarrier[3];
        QComb = 0.0;
        
    }
    else
    {
        QComb = -(dHeatSRCarrier[0]+dHeatSRCarrier[1]+dHeatSRCarrier[2]+dHeatSRCarrier[3])/dt;
    }
// dhsTrans = 0.0;
// QComb = 0.0;
//     const scalar xsi = min(T/cloud.constProps().TMax(), 1.0);
//     const scalar coeff =
//         (1.0 - xsi*xsi)*cloud.constProps().hRetentionCoeff();
// 
//     Sh += coeff*hReaction/dt;
// 
//     dhsTrans += (1.0 - coeff)*hReaction;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const ReactingMultiphaseParcel<ParcelType>& p
)
:
    ParcelType(p),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const ReactingMultiphaseParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelIO.C"

// ************************************************************************* //
