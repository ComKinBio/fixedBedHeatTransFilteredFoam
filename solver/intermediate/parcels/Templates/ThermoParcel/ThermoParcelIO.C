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

#include "ThermoParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ThermoParcel<ParcelType>::propertyList_ =
    Foam::ThermoParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::ThermoParcel<ParcelType>::sizeofFields_
(
    sizeof(ThermoParcel<ParcelType>)
  - offsetof(ThermoParcel<ParcelType>, T_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    T_(0.0),
    Cp_(0.0),
    kp_(0.0),
    dT_(0.0),    
    ppT4sum_(0.0),
    neighbourNum_(0.0),
    neighbourMax_(0.0),
    neighbourMin_(0.0),
    dBeforeDEM_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            T_ = readScalar(is);
            Cp_ = readScalar(is);
            kp_ = readScalar(is);
            dT_ = readScalar(is);            
            ppT4sum_ = readScalar(is);
            neighbourNum_ = readScalar(is);
            neighbourMax_ = readScalar(is);
            neighbourMin_ = readScalar(is);
            dBeforeDEM_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&T_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "ThermoParcel::ThermoParcel(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::ThermoParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, T);

    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Cp);
    
    IOField<scalar> kp(c.fieldIOobject("kp", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, kp);
    
    IOField<scalar> dT(c.fieldIOobject("dT", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, dT);
    
    IOField<scalar> ppT4sum(c.fieldIOobject("ppT4sum", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, ppT4sum);
    
    IOField<scalar> neighbourNum(c.fieldIOobject("neighbourNum", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, neighbourNum);
    
    IOField<scalar> neighbourMax(c.fieldIOobject("neighbourMax", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, neighbourMax);
    
    IOField<scalar> neighbourMin(c.fieldIOobject("neighbourMin", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, neighbourMin);
    
    IOField<scalar> dBeforeDEM(c.fieldIOobject("dBeforeDEM", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, dBeforeDEM);


    label i = 0;
    forAllIter(typename Cloud<ThermoParcel<ParcelType>>, c, iter)
    {
        ThermoParcel<ParcelType>& p = iter();

        p.T_ = T[i];
        p.Cp_ = Cp[i];
        p.kp_ = kp[i];
        p.dT_ = dT[i];
        p.ppT4sum_ = ppT4sum[i];
        p.neighbourNum_ = neighbourNum[i];
        p.neighbourMax_ = neighbourMax[i];
        p.neighbourMin_ = neighbourMin[i];
        p.dBeforeDEM_ = dBeforeDEM[i];
    
        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ThermoParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::NO_READ), np);
    IOField<scalar> kp(c.fieldIOobject("kp", IOobject::NO_READ), np);
    IOField<scalar> dT(c.fieldIOobject("dT", IOobject::NO_READ), np);
    IOField<scalar> ppT4sum(c.fieldIOobject("ppT4sum", IOobject::NO_READ), np);
    IOField<scalar> neighbourNum(c.fieldIOobject("neighbourNum", IOobject::NO_READ), np);
    IOField<scalar> neighbourMax(c.fieldIOobject("neighbourMax", IOobject::NO_READ), np);
    IOField<scalar> neighbourMin(c.fieldIOobject("neighbourMin", IOobject::NO_READ), np);
    IOField<scalar> dBeforeDEM(c.fieldIOobject("dBeforeDEM", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<ThermoParcel<ParcelType>>, c, iter)
    {
        const ThermoParcel<ParcelType>& p = iter();

        T[i] = p.T_;
        Cp[i] = p.Cp_;
        kp[i] = p.kp_;
        dT[i] = p.dT_;
        ppT4sum[i] = p.ppT4sum_;
        neighbourNum[i] = p.neighbourNum_;
        neighbourMax[i] = p.neighbourMax_;
        neighbourMin[i] = p.neighbourMin_;
        dBeforeDEM[i] = p.dBeforeDEM_;
        i++;
    }

    T.write(np > 0);
    Cp.write(np > 0);
    kp.write(np > 0);
    dT.write(np > 0);
    ppT4sum.write(np > 0);
    neighbourNum.write(np > 0);
    neighbourMax.write(np > 0);
    neighbourMin.write(np > 0);
    dBeforeDEM.write(np > 0);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ThermoParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.T()
            << token::SPACE << p.Cp()
            << token::SPACE << p.kp()
            << token::SPACE << p.dT()
            << token::SPACE << p.ppT4sum()
            << token::SPACE << p.neighbourNum()
            << token::SPACE << p.neighbourMax()
            << token::SPACE << p.neighbourMin()
            << token::SPACE << p.dBeforeDEM();

    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.T_),
            ThermoParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const ThermoParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
