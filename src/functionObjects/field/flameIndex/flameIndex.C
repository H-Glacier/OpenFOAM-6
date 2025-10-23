/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025
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
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

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
    if
    (
        fieldNames_.size() >= 2
     && foundObject<volScalarField>(fieldNames_[0])
     && foundObject<volScalarField>(fieldNames_[1])
    )
    {
        const volScalarField& YF = lookupObject<volScalarField>(fieldNames_[0]);
        const volScalarField& YO = lookupObject<volScalarField>(fieldNames_[1]);

        const tmp<volVectorField> tGradYF = fvc::grad(YF);
        const tmp<volVectorField> tGradYO = fvc::grad(YO);
        const volVectorField& gradYF = tGradYF();
        const volVectorField& gradYO = tGradYO();

        const tmp<volScalarField> tmagGradYF = mag(gradYF);
        const tmp<volScalarField> tmagGradYO = mag(gradYO);
        const volScalarField& magGradYF = tmagGradYF();
        const volScalarField& magGradYO = tmagGradYO();

        // Avoid division by zero using SMALL
        tmp<volScalarField> tFI
        (
            (gradYF & gradYO) / (max(magGradYF, SMALL) * max(magGradYO, SMALL))
        );

        return store(resultName_, tFI);
    }

    return false;
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
    setResultName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::flameIndex::~flameIndex()
{}


// ************************************************************************* //
