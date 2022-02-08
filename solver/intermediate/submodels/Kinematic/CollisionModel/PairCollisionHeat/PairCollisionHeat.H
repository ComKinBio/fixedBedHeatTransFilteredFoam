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
    Foam::PairCollisionHeat

Description

SourceFiles
    PairCollisionHeat.C

\*---------------------------------------------------------------------------*/

#ifndef PairCollisionHeat_H
#define PairCollisionHeat_H

#include "CollisionModel.H"
#include "InteractionLists.H"
#include "WallSiteData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class PairModel;

template<class CloudType>
class WallModel;


/*---------------------------------------------------------------------------*\
                     Class PairCollisionHeat Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PairCollisionHeat
:
    public CollisionModel<CloudType>
{
    // Static data

        //- Tolerance to determine flat wall interactions
        static scalar cosPhiMinFlatWall;

        //- Distance to multiply the distance to a wall with within
        //  which a duplicate flat wall interaction is defined.  If
        //  cosPhiFlatWall = 1 - e, this = sqrt(2e - e^2) < sqrt(3e)
        //  for all e < 1, giving a safe tolerance (i.e. no double
        //  interactions next to each other.)
        static scalar flatWallDuplicateExclusion;


    // Private Data

        //- PairModel to calculate the interaction between two parcels
        autoPtr<PairModel<CloudType>> pairModel_;

        //- WallModel to calculate the interaction between the parcel and walls
        autoPtr<WallModel<CloudType>> wallModel_;

        //- Interactions lists determining which cells are in
        //  interaction range of each other
        InteractionLists<typename CloudType::parcelType> il_;


    // Private Member Functions

        //- Pre collision tasks
        void preInteraction();
        
        //- Pre collision tasks for particle-particle heat transfer
        void preInteractionHeat();

        //- Interactions between parcels
        void parcelInteraction();
        
        //- Interactions between parcels for particle-particle heat transfer
        void parcelInteractionHeat();

        //- Interactions between real (on-processor) particles
        void realRealInteraction();

        //- Interactions between real and referred (off processor) particles
        void realReferredInteraction();
        
        //- Interactions between real (on-processor) particles for particle-particle heat transfer
        void realRealInteractionHeat();

        //- Interactions between real and referred (off processor) particles for particle-particle heat transfer
        void realReferredInteractionHeat();

        //- Interactions with walls
        void wallInteraction();
        
        //- Interactions with walls heat transfer
        void wallInteractionHeat();

        bool duplicatePointInList
        (
            const DynamicList<point>& existingPoints,
            const point& pointToTest,
            scalar duplicateRangeSqr
        ) const;

        bool duplicatePointInList
        (
            const DynamicList<point>& existingPoints,
            const point& pointToTest,
            const scalarList& duplicateRangeSqr
        ) const;

        //- Post collision tasks
        void postInteraction();

        //- Calculate the pair force between parcels
        void evaluatePair
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const;
        
        //- Calculate the pair heat transfer between parcels
        void evaluatePairHeat
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const;

        //- Calculate the wall forces on a parcel
        void evaluateWall
        (
            typename CloudType::parcelType& p,
            const List<point>& flatSitePoints,
            const List<WallSiteData<vector>>& flatSiteData,
            const List<point>& sharpSitePoints,
            const List<WallSiteData<vector>>& sharpSiteData
        ) const;
        
        //- Calculate  heat transfer between parcels and wall
        void evaluateWallHeat
        (
            typename CloudType::parcelType& p,
            const List<point>& flatSitePoints,
            const List<WallSiteData<vector>>& flatSiteData,
            const List<point>& sharpSitePoints,
            const List<WallSiteData<vector>>& sharpSiteData
        ) const;


public:

    //- Runtime type information
    TypeName("pairCollisionHeat");


    // Constructors

        //- Construct from components
        PairCollisionHeat(const dictionary& dict, CloudType& owner);

        //- Construct copy
        PairCollisionHeat(const PairCollisionHeat<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CollisionModel<CloudType>> clone() const
        {
            return autoPtr<CollisionModel<CloudType>>
            (
                new PairCollisionHeat<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PairCollisionHeat();


    // Member Functions

        //- Return the number of times to subcycle the current
        //  timestep to meet the criteria of the collision model.
        virtual label nSubCycles() const;

        // Collision function
        virtual void collide();
        
        // Collision function for particle=particle heat transfer
        virtual void collideHeat();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PairCollisionHeat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
