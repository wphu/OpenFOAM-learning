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
    Foam::PICParcel

Description
    PIC parcel class

SourceFiles
    PICParcelI.H
    PICParcel.C
    PICParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef PICParcel_H
#define PICParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"
#include "PICCloud.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class PICParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const PICParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class PICParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class PICParcel
:
    public ParcelType
{
    // Private member data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold PIC particle constant properties
    class constantProperties
    {
        // Private Data

            //- Particle mass [kg] (constant)
            scalar mass_;

            //- Particle hard sphere diameter [m] (constant)
            scalar d_;

            //- Internal degrees of freedom
            direction internalDegreesOfFreedom_;

            //- Viscosity index
            scalar omega_;

            //- Electrical charge
            scalar charge_;

    public:

        // Constructors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);


        // Member Functions

            //- Return const access to the particle mass [kg]
            inline scalar mass() const;

            //- Return const access to the hard sphere diameter [m]
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the internalDegreesOfFreedom
            inline direction internalDegreesOfFreedom() const;

            //- Return the viscosity index
            inline scalar omega() const;

            //- Return the electrical charge
            inline scalar charge() const;

    };


    //- Use base tracking data
    //typedef typename ParcelType::trackingData trackingData;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public ParcelType::trackingData
    {
        // Interpolators for continuous phase fields


            const interpolationCellPoint<vector>& EInterp_;
            const interpolationCellPoint<vector>& BInterp_;

    public:

        // Constructors

            inline trackingData
            (
                const PICCloud<PICParcel<ParcelType>>& pc,
                const interpolationCellPoint<vector>& EInterp,
                const interpolationCellPoint<vector>& BInterp
            );

        // Member Functions

            inline const interpolationCellPoint<vector>& EInterp() const;
            inline const interpolationCellPoint<vector>& BInterp() const;
    };





protected:

    // Protected member data

        // Parcel properties

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Internal energy of the Parcel, covering all non-translational
            //  degrees of freedom [J]
            scalar Ei_;

            //- Parcel type id
            label typeId_;


public:

    //- Runtime type information
    TypeName("PICParcel");

    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from components
        inline PICParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const vector& U,
            const scalar Ei,
            const label typeId
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        inline PICParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const vector& U,
            const scalar Ei,
            const label typeId
        );

        //- Construct from Istream
        PICParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new PICParcel<ParcelType>(*this));
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

            autoPtr<PICParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<PICParcel<ParcelType>>
                (
                    new PICParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return type id
            inline label typeId() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to internal energy
            inline scalar Ei() const;


        // Edit

            //- Return access to velocity
            inline vector& U();

            //- Return access to internal energy
            inline scalar& Ei();


        // Main calculation loop

            // Tracking

                //- update velocity
                template<class TrackCloudType>
                bool updateVelocity
                (
                    TrackCloudType& cloud,
                    trackingData& td,
                    const scalar trackTime
                );

                //- Move the parcel
                template<class TrackCloudType>
                bool move
                (
                    TrackCloudType& cloud,
                    trackingData& td,
                    const scalar trackTime
                );


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackCloudType>
            bool hitPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackCloudType>
            void hitProcessorPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackCloudType>
            void hitWallPatch(TrackCloudType&, trackingData&);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


        // I-O

            static void readFields(Cloud<PICParcel<ParcelType>>& c);

            static void writeFields(const Cloud<PICParcel<ParcelType>>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const PICParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PICParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PICParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
