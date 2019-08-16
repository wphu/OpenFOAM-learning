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

#include "PICParcel.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::PICParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = cloud.pMesh();


	
    const volVectorField& BField = mesh.lookupObject<volVectorField>("B");
    vector B = BField.internalField()[p.cell()];

    const volVectorField& EField = mesh.lookupObject<volVectorField>("E");
    vector E = EField.internalField()[p.cell()];

    const constantProperties& constProps(cloud.constProps(typeId_));

    scalar C = trackTime * constProps.charge() / constProps.mass();
    vector dv1 = C * ((U_ ^ B) + E);
    vector dv2 = C * (((U_ + dv1 / 2) ^ B) + E);
    vector dv3 = C * (((U_ + dv2 / 2) ^ B) + E);
    vector dv4 = C * (((U_ + dv3) ^ B) + E);
    U_ += (dv1 + 2*dv2 + 2*dv3 + dv4)/6;


    // For reduced-D cases, the velocity used to track needs to be
    // constrained, but the actual U_ of the parcel must not be
    // altered or used, as it is altered by patch interactions an
    // needs to retain its 3D value for collision purposes.
    vector Utracking = U_;

    while (td.keepParticle && !td.switchProcessor && p.stepFraction() < 1)
    {

        Utracking = U_;

        // Apply correction to velocity to constrain tracking for
        // reduced-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre();

        const scalar f = 1 - p.stepFraction();
        p.trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::PICParcel<ParcelType>::hitPatch(TrackCloudType&, trackingData&)
{
    return false;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::hitProcessorPatch
(
    TrackCloudType&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::hitWallPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    const label wppIndex = this->patch();

    const wallPolyPatch& wpp =
        static_cast<const wallPolyPatch&>
        (
            this->mesh().boundaryMesh()[wppIndex]
        );

    const label wppLocalFace = wpp.whichFace(this->face());

    const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

    const scalar deltaT = cloud.pMesh().time().deltaTValue();

    const constantProperties& constProps(cloud.constProps(typeId_));

    scalar m = constProps.mass();

    scalar charge = constProps.charge();

    // energy
    scalar E = 0.5*m*(U_ & U_);

    cloud.particleFluxBF(typeId_)[wppIndex][wppLocalFace] += cloud.nParticle()/(deltaT*fA);
    cloud.heatFluxBF(typeId_)[wppIndex][wppLocalFace] += cloud.nParticle()*E/(deltaT*fA);

    td.keepParticle = false;
}


template<class ParcelType>
void Foam::PICParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);
    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::PICParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "PICParcelIO.C"

// ************************************************************************* //
