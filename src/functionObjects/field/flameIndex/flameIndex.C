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

\*---------------------------------------------------------------------------*/

#include "flameIndex.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(flameIndex, 0);
    addToRunTimeSelectionTable(functionObject, flameIndex, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::flameIndex::calc()
{
    // Expect exactly two scalar fields (fuel and oxidiser mass fractions)
    if
    (
        !foundObject<volScalarField>(fieldNames_[0])
     || !foundObject<volScalarField>(fieldNames_[1])
    )
    {
        return false;
    }

    const volScalarField& YF = lookupObject<volScalarField>(fieldNames_[0]);
    const volScalarField& YO = lookupObject<volScalarField>(fieldNames_[1]);

    const volVectorField gYF = fvc::grad(YF);
    const volVectorField gYO = fvc::grad(YO);

    // Compute normalised dot product with stabilised denominators
    const volScalarField maggYF = mag(gYF);
    const volScalarField maggYO = mag(gYO);

    const dimensionedScalar smallG
    (
        "smallG", dimless, SMALL
    );

    tmp<volScalarField> tFI
    (
        new volScalarField
        (
            IOobject
            (
                resultName_,
                time_.timeName(),
                mesh_
            ),
            (gYF & gYO)
          / (stabilise(maggYF, smallG) * stabilise(maggYO, smallG))
        )
    );

    return store(resultName_, tFI);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::flameIndex::flameIndex
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldsExpression(name, runTime, dict)
{
    setResultName("flameIndex");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::flameIndex::~flameIndex()
{}


// ************************************************************************* //
