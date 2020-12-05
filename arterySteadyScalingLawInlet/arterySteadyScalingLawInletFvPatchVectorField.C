/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "arterySteadyScalingLawInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::arterySteadyScalingLawInletFvPatchVectorField::
arterySteadyScalingLawInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    areaTotal_(gSum(patch().magSf())),
    areaAssumedDiameter_(2*std::sqrt(areaTotal_/PI)),
    inletVelocity_(1.43*pow(areaAssumedDiameter_, 2.55)/areaTotal_)
{}


Foam::arterySteadyScalingLawInletFvPatchVectorField::
arterySteadyScalingLawInletFvPatchVectorField
(
    const arterySteadyScalingLawInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    areaTotal_(gSum(patch().magSf())),
    areaAssumedDiameter_(2*std::sqrt(areaTotal_/PI)),
    inletVelocity_(1.43*pow(areaAssumedDiameter_, 2.55)/areaTotal_)
{}


Foam::arterySteadyScalingLawInletFvPatchVectorField::
arterySteadyScalingLawInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    areaTotal_(gSum(patch().magSf())),
    areaAssumedDiameter_(2*std::sqrt(areaTotal_/PI)),
    inletVelocity_(1.43*pow(areaAssumedDiameter_, 2.55)/areaTotal_)
{
    // Note: No need to have an updateCoeff method as the inlet value is constant 
    fvPatchVectorField::operator=(-patch().nf()*inletVelocity_);
}


Foam::arterySteadyScalingLawInletFvPatchVectorField::
arterySteadyScalingLawInletFvPatchVectorField
(
    const arterySteadyScalingLawInletFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    areaTotal_(pivpvf.areaTotal_),
    areaAssumedDiameter_(pivpvf.areaAssumedDiameter_),
    inletVelocity_(pivpvf.inletVelocity_)
{}


Foam::arterySteadyScalingLawInletFvPatchVectorField::
arterySteadyScalingLawInletFvPatchVectorField
(
    const arterySteadyScalingLawInletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    areaTotal_(pivpvf.areaTotal_),
    areaAssumedDiameter_(pivpvf.areaAssumedDiameter_),
    inletVelocity_(pivpvf.inletVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// extract the normals and the total area of the inlet patch

void Foam::arterySteadyScalingLawInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::arterySteadyScalingLawInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::arterySteadyScalingLawInletFvPatchVectorField::write(Ostream& os) const
{
    // scalar areaTotal_;
    // // the area is in metres calculated by the constructors
    // scalar areaAssumedDiameter_;
    // // the flow rate is in m^3/s
    // scalar inletVelocity_;
    fvPatchVectorField::write(os);
    os.writeKeyword("area") << areaTotal_ << token::END_STATEMENT << nl;
    os.writeKeyword("areaAssumedDiameter") << areaAssumedDiameter_ << token::END_STATEMENT << nl;
    os.writeKeyword("velocity") << inletVelocity_ << token::END_STATEMENT << nl;
    // writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        arterySteadyScalingLawInletFvPatchVectorField
    );
}

// ************************************************************************* //
