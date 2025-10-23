/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionObjects::flameIndex

Description
    Calculates the Takeno flame index
        FI = (grad(Y_F) . grad(Y_O)) / (|grad(Y_F)| |grad(Y_O)|)

    The operation expects two scalar volume fields (fuel mass fraction and
    oxidiser mass fraction) supplied via the 'fields' entry in the function
    object dictionary, e.g.:

    
    flameIndex1
    {
        type    flameIndex;
        libs    ("libfieldFunctionObjects.so");
        fields  (CH4 O2);    // fuel, oxidiser
        result  FI;          // optional output name
        executeControl  writeTime;
        writeControl    writeTime;
    }
    

See also
    Foam::functionObjects::fieldsExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    flameIndex.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_flameIndex_H
#define functionObjects_flameIndex_H

#include "fieldsExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class flameIndex Declaration
\*---------------------------------------------------------------------------*/

class flameIndex
:
    public fieldsExpression
{
    // Private Member Functions

        //- Calculate the flame index field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("flameIndex");


    // Constructors

        //- Construct from Time and dictionary
        flameIndex
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~flameIndex();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
