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
    Foam::ReactingMultiphaseParcel

Description
    Multiphase variant of the reacting parcel class with one/two-way coupling
    with the continuous phase.

SourceFiles
    ReactingMultiphaseParcelI.H
    ReactingMultiphaseParcel.C
    ReactingMultiphaseParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseParcel_H
#define ReactingMultiphaseParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingMultiphaseParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingMultiphaseParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                 Class ReactingMultiphaseParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingMultiphaseParcel
:
    public ParcelType
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    // IDs of phases in ReacingParcel phase list (Y)

        static const label GAS;
        static const label LIQ;
        static const label SLD;


    //- Class to hold reacting multiphase particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private Data

            //- Devolatilisation activation temperature [K]
            demandDrivenEntry<scalar> TDevol_;

            //- Latent heat of devolatilisation [J/kg]
            demandDrivenEntry<scalar> LDevol_;

            //- Fraction of enthalpy retained by parcel due to surface
            //  reactions
            demandDrivenEntry<scalar> hRetentionCoeff_;

            //- coarse grid flag
            demandDrivenEntry<bool> coarseGrid_;
            
            //- max solid phase volume fraction
            demandDrivenEntry<scalar> alphaMax_;   
            
            //- disable htc flag
            demandDrivenEntry<bool> disableHtc_;
            
            //- disable combustion flag
            demandDrivenEntry<bool> noCombustion_;
            
            //- use combustion heat ratio between solid and gas flag
            demandDrivenEntry<bool> heatRatioFlag_;

            //- use heat conduction between particles flag
            demandDrivenEntry<bool> ppConductionFlag_;
            
            //- combustion control flag
            demandDrivenEntry<bool> controledCombustion_;
            
            //- combusiton occurs after the degree of Devo
            demandDrivenEntry<scalar> combustionAt_;
            
            //- surface temperature quation solver control flag
            demandDrivenEntry<bool> surfaceExplicit_;
            
            //- surface species concentrations sampling control flag
            demandDrivenEntry<bool> coarseGridSurfaceCombustion_;
            
            //- correct radiaon heat for Devo particle
            demandDrivenEntry<scalar> radiationCorrection_;
            
            //- enable conductive radiaon
            demandDrivenEntry<bool> ppRadiationFlag_;
            
            //- shrinkageMode selection, 1 is default one ash layer have constant porosity as 0.65;
            //- 0 is calculated porosity of ash layer assuming ash density is 2000kg/m-3
            demandDrivenEntry<label> shrinkageMode_;
            
            //- shrinkageMode factor for char -> ash
            demandDrivenEntry<scalar> shrinkageFactorAsh_;
            
// SC * * * * * * * * * * * * * * * Variables Declaration* * * * * * * * * * * * * * * * * * * * * * //           
            //- Layer initial properties
            demandDrivenEntry<scalar> Tb00_;
            demandDrivenEntry<scalar> Tb10_;
            demandDrivenEntry<scalar> Tb20_;
            demandDrivenEntry<scalar> Tb30_;
            demandDrivenEntry<scalar> rb00_;
            demandDrivenEntry<scalar> rb10_;
            demandDrivenEntry<scalar> rb20_;    
            demandDrivenEntry<scalar> rb30_;
            demandDrivenEntry<scalar> Tp00_;
            demandDrivenEntry<scalar> Tp10_;     
            demandDrivenEntry<scalar> Tp20_;  
            demandDrivenEntry<scalar> Tp30_;
            demandDrivenEntry<scalar> deltaTime_;
            demandDrivenEntry<scalar> hCoeff_;         
// -------------------------Variables in class constantProperties Public-----------------------------// 


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access

            //- Return const access to the devolatilisation temperature
            inline scalar TDevol() const;

            //- Return const access to the latent heat of devolatilisation
            inline scalar LDevol() const;

            //- Return const access to the fraction of enthalpy retained by
            //  parcel due to surface reactions
            inline scalar hRetentionCoeff() const;
            
            //- Return const access to the coarse grid flag
            inline bool coarseGrid() const;
            
            //- Return const access to the max solid phase volume fraction
            inline scalar alphaMax() const;
            
            //- Return const access to the diable htc flag
            inline bool disableHtc() const;
            
            //- Return const access to the diable htc flag
            inline bool noCombustion() const;
            
            //- Return const access to the diable htc flag
            inline bool heatRatioFlag() const;
            
            //- Return const access to heat conduction between particles flag
            inline bool ppConductionFlag() const;
            
            //- Return const access to combustion control flag
            inline bool controledCombustion() const;
            
            //- Return const access to combusiton occurs after the degree of Devo
            inline scalar combustionAt() const;
            
            //- surface temperature quation solver control flag
            inline bool surfaceExplicit() const;
            
            //-  Return const access to surface species concentrations sampling control flag
            inline bool coarseGridSurfaceCombustion() const;
            
            //- Return const access to correct radiaon heat for Devo particle
            inline scalar radiationCorrection() const;
            
            //- enable conductive radiaon
            inline bool ppRadiationFlag() const;
            
            //- Return const access to shrinkageMode
            inline label shrinkageMode() const;
            
            //- Return const access to shrinkageMode factor for char -> ash
            inline scalar shrinkageFactorAsh() const;
            
// SC * * * * * * * * * * * * * * * Variables get function* * * * * * * * * * * * * * * * * * * * * * //  
            //- Return const access to layer inital properties
            inline scalar Tb00() const;
            inline scalar Tb10() const;
            inline scalar Tb20() const;
            inline scalar Tb30() const;
            inline scalar rb00() const;
            inline scalar rb10() const;
            inline scalar rb20() const;
            inline scalar rb30() const;
            inline scalar Tp00() const;
            inline scalar Tp10() const;
            inline scalar Tp20() const;
            inline scalar Tp30() const;

            //particle layer delta time
            inline scalar deltaTime() const;

            //heat transfer coefficient
            inline scalar hCoeff() const;          
// -------------------------Variables in class constantProperties Public-----------------------------//     
	};


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;


private:

    // Private Member Functions

        //- Return the mixture effective specific heat capacity
        template<class TrackCloudType>
        scalar CpEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective sensible enthalpy
        template<class TrackCloudType>
        scalar HsEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective latent heat
        template<class TrackCloudType>
        scalar LEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Update the mass fractions (Y, YGas, YLiquid, YSolid)
        scalar updateMassFractions
        (
            const scalar mass0,
            const scalarField& dMassGas,
            const scalarField& dMassLiquid,
            const scalarField& dMassSolid
        );
        

// SC * * * * * * * * * * * * * * * sub-model function* * * * * * * * * * * * * * * * * * * * * * //         
     // Particle layer Functions
        
        //- get radius of circle surface from the volume of cylinderL
        scalar radiusForCylinder
        (
            const scalar Xi, 
            const scalar V0
        );
      
        //- Particle surface area sphere
        scalar Area_Sph
        (
            const scalar radius
        );
        
        
        //- Particle surface area cylinderL
        scalar Area_cylinderL
        (
            const scalar radius,
            const scalar Xi
        );
        
        //- Particle volume cylinderL
        scalar Volume_cylinderL
        (
            const scalar radius,
            const scalar Xi
        );
        
        //- Shell volume sphere
        scalar Vol_Rin 
        (
            const scalar rin,
            const scalar rout
        );
        
        //- Shell volume cylinderL
        scalar Vol_Rin_cylinderL
        (
            const scalar rin,
            const scalar rout,
            const scalar Xi
        );
        
        //- Radius of the shell sphere
        scalar R_Par
        (
            const scalar rin,
            const scalar rout
        );
        
        //- Radius of the shell cylinderL
        scalar R_Par_cylinderL
        (
            const scalar rin,
            const scalar rout,
            const scalar Xi
        );
        
        //- Compute A*kp*d/dr sphere
        scalar d_dr
        (
            const scalar kp,
            const scalar Ab,
            const scalar ri,
            const scalar rj
        );
        
        //- Compute A*kp*d/dr cylinderL
        scalar d_dr_cylinderL
        (
            const scalar kp,
            const scalar Ab,
            const scalar ri,
            const scalar rj,
            const scalar Xi
        );
        
        //- Function Fb1 based on Antoine expression 
        scalar Fb1 
        (
            const scalar Tb1
        );
        
        //- Heat balance over the outer particle surface
        bool eq4_check
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
        );
        
        //- eq4 Henrik C&F for cylinderL
        bool eq4_cylinderL_check
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
        );
        
        
        //- eq4 Henrik C&F
        scalar eq4
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
        );
        
        //- eq4 Henrik C&F
        scalar eq4_Explicit
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
        );
        
         //- Heat balance over the boundary of the layer
         //- Derived from eq7 Henrik C&F
        scalar eq7_2
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
        );
        
        scalar eq7_3
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
        );
        
        //- eq4 Henrik C&F for cylinderL
        scalar eq4_cylinderL
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
        );
        
        //- eq4 Henrik C&F for cylinderL
        scalar eq4_cylinderL_Explicit
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
        );
        
         //- Heat balance over the boundary of the layer
         //- Derived from eq7 Henrik C&F for cylinderL
        scalar eq7_2_cylinderL
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
        );
        
        scalar eq7_3_cylinderL
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
        );
        //- Heat capacity of the soild
        
        scalar cp_p
        (
            const label layer,
            const scalar Tp,
            const scalar moist_WB
        );
        
        scalar cp_p_modified
        (
            const label layer,
            const scalar Tp,
            const scalar moist_WB
        );
        
        //- Density of the soild
        
        scalar rho_p
        (
            const label layer,
            const scalar Tp,
            const scalar moist_WB
        );
        
        //- Conductivity of the soild
        
        scalar kp_p
        (
            const label layer,
            const scalar Tp
        );
        
        //- Heat capacity of the water vapor
        scalar cp_water_vapor
        (
            const scalar T
        );
        scalar deltaHvap
        (
            const scalar T
        );
// -------------------------sub-model function-----------------------------//        


protected:

    // Protected data 
    


// SC * * * * * * * * * * * * * * * Variables during SC * * * * * * * * * * * * * * * * * * * * * * //        
        // Boundary value
            
            //- Tb_ tempearture [K], rb_ radius [m]
            scalar Tb0_, Tb1_, Tb2_, Tb3_, rb0_, rb1_, rb2_, rb3_;
            
            //- Particle value
            // Tp_ tempearture [K], rp_ mass-weighted radius [m], layer mass [k]
            scalar Tp0_, Tp1_, Tp2_, Tp3_, mp0_, mp1_, mp2_, mp3_;
            
            //- Mass loss/gain due to drying and devolatilization and char combustion (kg)
            scalar rDry_, rDevo_, rChar_, rComb_, ash_inchar_t_;
            
            //- Heat flux due to drying and combustion [W]
            scalar QDry_, QComb_;
            
            //- Flag for boiling and devolatilisation
            
            label flagBoiling_, flagDevo_;
            
            //- Accumulated time [s]
            scalar cumTime_;
            
            //- Qin, QLayer 
            scalar Qin_, QLayer_;

            //variables in debug 
            scalar VID1_, VID2_, VID3_, VID4_;
// -------------------------Variables during SC in parcel class protected-----------------------------//  


        // Parcel properties

            //- Mass fractions of gases []
            scalarField YGas_;

            //- Mass fractions of liquids []
            scalarField YLiquid_;

            //- Mass fractions of solids []
            scalarField YSolid_;

            //- Flag to identify if the particle can devolatilise and combust
            //  Combustion possible only after volatile content falls below
            //  threshold value.  States include:
            //  0 = can devolatilise, cannot combust but can change
            //  1 = can devolatilise, can combust
            // -1 = cannot devolatilise or combust, and cannot change
            label canCombust_;


    // Protected Member Functions

        //- Calculate Devolatilisation
        template<class TrackCloudType>
        void calcDevolatilisation
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar age,          // age
            const scalar Ts,           // surface temperature
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const scalar mass0,        // mass (initial on injection)
            const scalarField& YGasEff,// gas component mass fractions
            const scalarField& YLiquidEff,// liquid component mass fractions
            const scalarField& YSolidEff,// solid component mass fractions
            label& canCombust,          // 'can combust' flag
            scalarField& dMassDV,      // mass transfer - local to particle
//SC * * * * * * * * * * * * * * * add solid mass due to Devo models ------------------------- //  
            scalarField& dMassSOLID,   // soild mass update         
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        ) const;
        
        //- Calculate surface reactions
        template<class TrackCloudType>
        void calcSurfaceReactions
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar d,            // radius
            const scalar di,            // inner layer radius
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label canCombust,     // 'can combust' flag
            const scalar N,            // flux of species emitted from particle
            const scalarField& YMix,   // mixture mass fractions
            const scalarField& YGas,   // gas-phase mass fractions
            const scalarField& YLiquid,// liquid-phase mass fractions
            const scalarField& YSolid, // solid-phase mass fractions
            scalarField& dMassSRGas,   // gas-phase mass transfer - local
            scalarField& dMassSRLiquid,// liquid-phase mass transfer - local
            scalarField& dMassSRSolid, // solid-phase mass transfer - local
            scalarField& dMassSRCarrier, // carrier phase mass transfer
            scalar& Sh,                // explicit particle enthalpy source
            scalar& dhsTrans,           // sensible enthalpy transfer to carrier
            scalar& QComb,               // record combustion heat
            const scalar Re,            // Re number
            const scalar Tc,            // continue phase temperature
            const scalar rhoc,          // continue phase density
            const scalar muc,              // continue phase viscosity
            const scalar Xi,            //cylinder shape const
            scalarField& dHeatSRCarrier,   //SR heat
            const label particleShape,
            const scalarField& Csavg,     //carriers surface concentrations and bed voidage
            const bool coarseGrid,
            const scalar deq,           //cylinder eq d
            const bool heatRatio,
            const bool limitedCombustion
        ) const;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("ReactingMultiphaseParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " nGas(Y1..YN)"
          + " nLiquid(Y1..YN)"
          + " nSolid(Y1..YN)"
        );


    // Constructors

        //- Construct from mesh, position and topology
        //  Other properties initialised as null
        inline ReactingMultiphaseParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ReactingMultiphaseParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline ReactingMultiphaseParcel
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
        );

        //- Construct from Istream
        ReactingMultiphaseParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ReactingMultiphaseParcel(const ReactingMultiphaseParcel& p);

        //- Construct as a copy
        ReactingMultiphaseParcel
        (
            const ReactingMultiphaseParcel& p,
            const polyMesh& mesh
        );

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ReactingMultiphaseParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ReactingMultiphaseParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ReactingMultiphaseParcel<ParcelType>> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<ReactingMultiphaseParcel<ParcelType>>
                (
                    new ReactingMultiphaseParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to mass fractions of gases
            inline const scalarField& YGas() const;

            //- Return const access to mass fractions of liquids
            inline const scalarField& YLiquid() const;

            //- Return const access to mass fractions of solids
            inline const scalarField& YSolid() const;

            //- Return const access to the canCombust flag
            inline label canCombust() const;
            
            


// SC * * * * * * * * * * * * * * * Variables access functions * * * * * * * * * * * * * * * * * * * * * * //        
            inline scalar Tb0() const;
            inline scalar Tb1() const;
            inline scalar Tb2() const;
            inline scalar Tb3() const;
            inline scalar rb0() const;
            inline scalar rb1() const;
            inline scalar rb2() const;
            inline scalar rb3() const;
            inline scalar Tp0() const;
            inline scalar Tp1() const;
            inline scalar Tp2() const;
            inline scalar Tp3() const;
            inline scalar mp0() const;
            inline scalar mp1() const;
            inline scalar mp2() const;
            inline scalar mp3() const;
            inline scalar rDry() const;
            inline scalar rDevo() const;
            inline scalar rChar() const;
            inline scalar rComb() const;
            inline scalar ash_inchar_t() const;
            inline scalar QDry() const; 
            inline scalar QComb() const; 
            inline label flagBoiling() const;
            inline label flagDevo() const;
            inline scalar cumTime() const;
            inline scalar Qin() const; 
            inline scalar QLayer() const; 
            inline scalar VID1() const; 
            inline scalar VID2() const; 
            inline scalar VID3() const; 
            inline scalar VID4() const; 
// -------------------------Variables access functions parcel class Public-----------------------------//  


        // Edit

            //- Return access to mass fractions of gases
            inline scalarField& YGas();

            //- Return access to mass fractions of liquids
            inline scalarField& YLiquid();

            //- Return access to mass fractions of solids
            inline scalarField& YSolid();

            //- Return access to the canCombust flag
            inline label& canCombust();


// SC * * * * * * * * * * * * * * * Variables edit functions * * * * * * * * * * * * * * * * * * * * * * // 
            inline scalar& Tb0();
            inline scalar& Tb1();
            inline scalar& Tb2();
            inline scalar& Tb3();
            inline scalar& rb0();
            inline scalar& rb1();
            inline scalar& rb2();
            inline scalar& rb3();
            inline scalar& Tp0();
            inline scalar& Tp1();
            inline scalar& Tp2();
            inline scalar& Tp3();
            inline scalar& mp0();
            inline scalar& mp1();
            inline scalar& mp2();
            inline scalar& mp3();
            inline scalar& rDry();
            inline scalar& rDevo();
            inline scalar& rChar();
            inline scalar& rComb();
            inline scalar& ash_inchar_t();
            inline scalar& QDry();
            inline scalar& QComb();
            inline label& flagBoiling();
            inline label& flagDevo();
            inline scalar& cumTime();
            inline scalar& Qin(); 
            inline scalar& QLayer(); 
            inline scalar& VID1(); 
            inline scalar& VID2(); 
            inline scalar& VID3(); 
            inline scalar& VID4(); 
// -------------------------Variables edit functions parcel class Public-----------------------------//  

            //- Return const access to the particle ash layer heat conductivity
            scalar kpp
            (
                const scalar Tp
            );

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );
            
            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calcV
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Read - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingMultiphaseParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingMultiphaseParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
