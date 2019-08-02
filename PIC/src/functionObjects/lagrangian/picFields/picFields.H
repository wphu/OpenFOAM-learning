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
    Foam::functionObjects::picFields

Description
    Calculate intensive fields:
    - UMean
    - translationalT
    - internalT
    - overallT
    from averaged extensive fields from a PIC calculation.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    picFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_picFields_H
#define functionObjects_picFields_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class picFields Declaration
\*---------------------------------------------------------------------------*/

class picFields
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("picFields");


    // Constructors

        //- Construct from Time and dictionary
        picFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        picFields(const picFields&) = delete;


    //- Destructor
    virtual ~picFields();


    // Member Functions

        //- Read the picFields data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Calculate and write the PIC fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const picFields&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
