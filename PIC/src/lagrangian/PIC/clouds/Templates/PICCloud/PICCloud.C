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

#include "PICCloud.H"
#include "BinaryCollisionModel.H"
#include "WallInteractionModel.H"
#include "InflowBoundaryModel.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"
#include "polyMeshTetDecomposition.H"
#include "AveragingMethod.H"
#include "interpolationCellPoint.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParcelType>
void Foam::PICCloud<ParcelType>::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< "    " << id << endl;

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] =
        typename ParcelType::constantProperties(molDict);
    }
}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(typename PICCloud<ParcelType>, *this, iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }
}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::initialise
(
    const IOdictionary& picInitialiseDict
)
{
    Info<< nl << "Initialising particles" << endl;

    const scalar temperature
    (
        readScalar(picInitialiseDict.lookup("temperature"))
    );

    const vector velocity(picInitialiseDict.lookup("velocity"));

    const dictionary& numberDensitiesDict
    (
        picInitialiseDict.subDict("numberDensities")
    );

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar> numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = readScalar
        (
            numberDensitiesDict.lookup(molecules[i])
        );
    }

    numberDensities /= nParticle_;

    forAll(mesh_.cells(), celli)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
            mesh_,
            celli
        );

        forAll(cellTets, tetI)
        {
            const tetIndices& cellTetIs = cellTets[tetI];
            tetPointRef tet = cellTetIs.tet(mesh_);
            scalar tetVolume = tet.mag();

            forAll(molecules, i)
            {
                const word& moleculeName(molecules[i]);

                label typeId(findIndex(typeIdList_, moleculeName));

                if (typeId == -1)
                {
                    FatalErrorInFunction
                        << "typeId " << moleculeName << "not defined." << nl
                        << abort(FatalError);
                }

                const typename ParcelType::constantProperties& cP =
                constProps(typeId);

                scalar numberDensity = numberDensities[i];

                // Calculate the number of particles required
                scalar particlesRequired = numberDensity*tetVolume;

                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > rndGen_.scalar01()
                )
                {
                    nParticlesToInsert++;
                }

                
                for (label pI = 0; pI < nParticlesToInsert; pI++)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U = equipartitionLinearVelocity
                    (
                        temperature,
                        cP.mass()
                    );

                    scalar Ei = equipartitionInternalEnergy
                    (
                        temperature,
                        cP.internalDegreesOfFreedom()
                    );

                    U += velocity;

                    addNewParcel(p, celli, U, Ei, typeId);
                }
            }
        }
    }


    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const typename ParcelType::constantProperties& cP = constProps
    (
        mostAbundantType
    );

    sigmaTcRMax_.primitiveFieldRef() = cP.sigmaT()*maxwellianMostProbableSpeed
    (
        temperature,
        cP.mass()
    );

    sigmaTcRMax_.correctBoundaryConditions();


}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::collisions()
{
    if (!binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    List<DynamicList<label>> subCells(8);

    scalar deltaT = mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;

    forAll(cellOccupancy_, celli)
    {
        const DynamicList<ParcelType*>& cellParcels(cellOccupancy_[celli]);

        label nC(cellParcels.size());

        if (nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            const point& cC = mesh_.cellCentres()[celli];

            forAll(cellParcels, i)
            {
                const ParcelType& p = *cellParcels[i];
                vector relPos = p.position() - cC;

                label subCell =
                    pos0(relPos.x()) + 2*pos0(relPos.y()) + 4*pos0(relPos.z());

                subCells[subCell].append(i);
                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            scalar sigmaTcRMax = sigmaTcRMax_[celli];

            scalar selectedPairs =
                collisionSelectionRemainder_[celli]
              + 0.5*nC*(nC - 1)*nParticle_*sigmaTcRMax*deltaT
               /mesh_.cellVolumes()[celli];

            label nCandidates(selectedPairs);
            collisionSelectionRemainder_[celli] = selectedPairs - nCandidates;
            collisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                label candidateP = rndGen_.sampleAB<label>(0, nC);

                // Declare the second collision candidate
                label candidateQ = -1;

                List<label> subCellPs = subCells[whichSubCell[candidateP]];
                label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        candidateQ = subCellPs[rndGen_.sampleAB<label>(0, nSC)];
                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        candidateQ = rndGen_.sampleAB<label>(0, nC);
                    } while (candidateP == candidateQ);
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = rndGen_.sampleAB<label>(0, nC);

                // // Select a possible second collision candidate
                // label candidateQ = rndGen_.sampleAB<label>(0, nC);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = rndGen_.sampleAB<label>(0, nC);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                ParcelType& parcelP = *cellParcels[candidateP];
                ParcelType& parcelQ = *cellParcels[candidateQ];

                scalar sigmaTcR = binaryCollision().sigmaTcR
                (
                    parcelP,
                    parcelQ
                );

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria because
                // the number of collision candidates selected was based on this

                if (sigmaTcR > sigmaTcRMax_[celli])
                {
                    sigmaTcRMax_[celli] = sigmaTcR;
                }

                if ((sigmaTcR/sigmaTcRMax) > rndGen_.scalar01())
                {
                    binaryCollision().collide
                    (
                        parcelP,
                        parcelQ
                    );

                    collisions++;
                }
            }
        }
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    sigmaTcRMax_.correctBoundaryConditions();

    if (collisionCandidates)
    {
        Info<< "    Collisions                      = "
            << collisions << nl
            << "    Acceptance rate                 = "
            << scalar(collisions)/scalar(collisionCandidates) << nl
            << endl;
    }
    else
    {
        Info<< "    No collisions" << endl;
    }
}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::resetFields()
{

    rhoN_ = dimensionedScalar( dimensionSet(0, -3, 0, 0, 0), vSmall);
    rhoQ_ =  dimensionedScalar( dimensionSet(0, -3, 1, 0, 0, 1, 0), vSmall);
    rhoQNew_ =  dimensionedScalar( dimensionSet(0, -3, 1, 0, 0, 1, 0), vSmall);

    forAll(rhoNSpeciesList_, iRhoNSpecies)
    {
        (*rhoNSpeciesList_[iRhoNSpecies]) = dimensionedScalar( dimensionSet(0, -3, 0, 0, 0), vSmall);
    }

    forAll(particleFluxList_, iParticleFlux)
    {
        (*rhoNSpeciesList_[iParticleFlux]) = dimensionedScalar( dimensionSet(0, -3, 0, 0, 0), vSmall);
    }

    forAll(heatFluxList_, iHeatFlux)
    {
        (*rhoNSpeciesList_[iHeatFlux]) = dimensionedScalar( dimensionSet(0, -3, 0, 0, 0), vSmall);
    }

    (*rhoQNewAverage_) = 0.0;
}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::calculateFields()
{
    scalarField& rhoN = rhoN_.primitiveFieldRef();
    scalarField& rhoQ = rhoQ_.primitiveFieldRef(); 
    scalarField& rhoQNew = rhoQNew_.primitiveFieldRef(); 

    forAllConstIter(typename PICCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label celli = p.cell();
        const label typeId = p.typeId();
        const scalar q = constProps(p.typeId()).charge();

        rhoN[celli]++;
        rhoQ[celli] += q;

        rhoNSpeciesList_[typeId]->primitiveFieldRef()[celli]++;

        //const tetIndices tetIs = p.currentTetIndices();
        //rhoQNewAverage_->add(p.coordinates(), tetIs, q);
    }

    rhoN *= nParticle_/mesh().cellVolumes();
    rhoN_.correctBoundaryConditions();

    rhoQ *= nParticle_/mesh().cellVolumes();
    rhoQ_.correctBoundaryConditions();

    forAll(rhoNSpeciesList_, iRhoNSpecies)
    {
        scalarField& rhoNSpecies = rhoNSpeciesList_[iRhoNSpecies]->primitiveFieldRef();
        rhoNSpecies *= nParticle_/mesh().cellVolumes();
        rhoNSpeciesList_[iRhoNSpecies]->correctBoundaryConditions();
    }

    rhoQNew = rhoQNewAverage_->averagePrimitiveField();
    rhoQNew *= nParticle_/mesh().cellVolumes();
    rhoQNew_.correctBoundaryConditions();


}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::PICCloud<ParcelType>::addNewParcel
(
    const vector& position,
    const label celli,
    const vector& U,
    const scalar Ei,
    const label typeId
)
{
    this->addParticle(new ParcelType(mesh_, position, celli, U, Ei, typeId));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::PICCloud<ParcelType>::PICCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    bool readFields
)
:
    Cloud<ParcelType>(mesh, cloudName, false),
    PICBaseCloud(),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(mesh_.nCells()),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    collisionSelectionRemainder_
    (
        IOobject
        (
            this->name() + ":collisionSelectionRemainder",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoQ_
    (
        IOobject
        (
            "rhoQ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoQNew_
    (
        IOobject
        (
            "rhoQNew",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    constProps_(),
    rndGen_(label(149382906) + 7183*Pstream::myProcNo()),
    boundaryT_
    (
        volScalarField
        (
            IOobject
            (
                "boundaryT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    boundaryU_
    (
        volVectorField
        (
            IOobject
            (
                "boundaryU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    binaryCollisionModel_
    (
        BinaryCollisionModel<PICCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        )
    ),
    wallInteractionModel_
    (
        WallInteractionModel<PICCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        )
    ),
    inflowBoundaryModel_
    (
        InflowBoundaryModel<PICCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        )
    )
{
    buildConstProps();
    buildCellOccupancy();


    rhoNSpeciesList_.setSize(typeIdList_.size());
    forAll(rhoNSpeciesList_, iRhoNSpecies)
    {
        word name ("rhoN_" + typeIdList_[iRhoNSpecies]);

        rhoNSpeciesList_[iRhoNSpecies].reset
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }


    particleFluxList_.setSize(typeIdList_.size());
    forAll(particleFluxList_, iParticleFlux)
    {
        word name ("particleFlux_" + typeIdList_[iParticleFlux]);

        particleFluxList_[iParticleFlux].reset
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    heatFluxList_.setSize(typeIdList_.size());
    forAll(heatFluxList_, iHeatFlux)
    {
        word name ("heatFlux_" + typeIdList_[iHeatFlux]);

        heatFluxList_[iHeatFlux].reset
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }


    
    dictionary AveragingMethod_dict
    (
        particleProperties_.subDict("moleculeProperties")
    );

    //typename ParcelType::trackingData td(*this);

    rhoQNewAverage_ = 
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                this->name() + ":rhoQNewAverage",
                this->db().time().timeName(),
                this->mesh()
            ),
            AveragingMethod_dict,
            this->mesh()
        )
    );
    
    



    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, i)
    {
        collisionSelectionRemainder_[i] = rndGen_.scalar01();
    }

    if (readFields)
    {
        ParcelType::readFields(*this);
    }
}


template<class ParcelType>
Foam::PICCloud<ParcelType>::PICCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& picInitialiseDict
)
    :
    Cloud<ParcelType>(mesh, cloudName, false),
    PICBaseCloud(),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( dimensionSet(0, 3, -1, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    collisionSelectionRemainder_
    (
        IOobject
        (
            this->name() + ":collisionSelectionRemainder",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    ),
    rhoN_
    (
        IOobject
        (
            this->name() + "rhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( dimensionSet(0, -3, 0, 0, 0), vSmall)
    ),
    rhoQ_
    (
        IOobject
        (
            this->name() + "rhoQ_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( dimensionSet(0, -3, 1, 0, 0, 1, 0), vSmall)
    ),
    rhoQNew_
    (
        IOobject
        (
            this->name() + "rhoQNew_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( dimensionSet(0, -3, 1, 0, 0, 1, 0), vSmall)
    ),
    constProps_(),
    rndGen_(label(971501) + 1526*Pstream::myProcNo()),
    boundaryT_
    (
        volScalarField
        (
            IOobject
            (
                "boundaryT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar( dimensionSet(0, 0, 0, 1, 0), 0)
        )
    ),
    boundaryU_
    (
        volVectorField
        (
            IOobject
            (
                "boundaryU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "zero",
                dimensionSet(0, 1, -1, 0, 0),
                Zero
            )
        )
    ),
    binaryCollisionModel_(),
    wallInteractionModel_(),
    inflowBoundaryModel_()
{
    clear();
    buildConstProps();
    initialise(picInitialiseDict);

    
    dictionary AveragingMethod_dict
    (
        particleProperties_.subDict("moleculeProperties")
    );

    //typename ParcelType::trackingData td(*this);

    rhoQNewAverage_ = 
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                this->name() + ":rhoQNewAverage",
                this->db().time().timeName(),
                this->mesh()
            ),
            AveragingMethod_dict,
            this->mesh()
        )
    );
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::PICCloud<ParcelType>::~PICCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::PICCloud<ParcelType>::evolve()
{
    //typename ParcelType::trackingData td(*this);

    // Reset the data collection fields
    resetFields();

    if (debug)
    {
        this->dumpParticlePositions();
    }

    // Insert new particles from the inflow boundary
    this->inflowBoundary().inflow();


    const volVectorField& E = mesh_.lookupObject<const volVectorField>("E");
    const volVectorField& B = mesh_.lookupObject<const volVectorField>("B");

    interpolationCellPoint<vector> EInterp(E);
    interpolationCellPoint<vector> BInterp(B);


    typename ParcelType::trackingData td(*this, EInterp, BInterp);

    // Move the particles ballistically with their current velocities
    Cloud<ParcelType>::move(*this, td, mesh_.time().deltaTValue());

    // Update cell occupancy
    buildCellOccupancy();

    // Calculate new velocities via stochastic collisions
    collisions();

    // Calculate the volume field data
    calculateFields();
}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::info() const
{
    label nPICParticles = this->size();
    reduce(nPICParticles, sumOp<label>());

    scalar nMol = nPICParticles*nParticle_;

    Info<< "Cloud name: " << this->name() << nl
        << "    Number of pic particles        = "
        << nPICParticles
        << endl;

    if (nPICParticles)
    {
        Info<< "    Number of molecules             = "
            << nMol << nl
            << endl;
    }
}


template<class ParcelType>
Foam::vector Foam::PICCloud<ParcelType>::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *rndGen_.sampleNormal<vector>();
}


template<class ParcelType>
Foam::scalar Foam::PICCloud<ParcelType>::equipartitionInternalEnergy
(
    scalar temperature,
    direction iDof
)
{
    // set to 0.0, or picFoam will encounter some weird errors
    // Maybe the method 'equipartitionInternalEnergy' and other relative codes should be removed in the future.
    return 0.0;
}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(typename PICCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        pObj<< "v " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << nl;
    }

    pObj.flush();
}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<ParcelType>::autoMap(mapper);

    // Update the cell occupancy field
    cellOccupancy_.setSize(mesh_.nCells());
    buildCellOccupancy();

    // Update the inflow BCs
    this->inflowBoundary().autoMap(mapper);
}


// ************************************************************************* //
