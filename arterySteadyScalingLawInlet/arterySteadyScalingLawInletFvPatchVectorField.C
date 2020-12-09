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
    area_(0.),
    areaAssumedDiameter_(0.),
    magnitudeVelocity_(0.),
    values_(p.size())
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
    area_(0.),
    areaAssumedDiameter_(0.),
    magnitudeVelocity_(0.),
    values_(p.size())
{

    if (dict.found("values")){
        dict.lookup("area") >> area_;
        dict.lookup("areaAssumedDiameter") >> areaAssumedDiameter_;
        dict.lookup("magnitudeVelocity") >> magnitudeVelocity_;
        vectorField data("values", dict, p.size());
        values_ = data;
        fvPatchVectorField::operator=(values_);
    } else {
        area_                = dict.lookupOrDefault("area", gSum(patch().magSf()));
        areaAssumedDiameter_ = dict.lookupOrDefault("areaAssumedDiameter",2*std::sqrt(area_/PI));
        const scalar flowRateTotal = 1.43*pow(areaAssumedDiameter_, 2.55);
        magnitudeVelocity_ = flowRateTotal/area_;
        values_             = -patch().nf()*magnitudeVelocity_;
        fvPatchVectorField::operator=(values_);
    }
}


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
    area_(ptf.area_),
    areaAssumedDiameter_(ptf.areaAssumedDiameter_),
    magnitudeVelocity_(ptf.magnitudeVelocity_),
    values_(ptf.values_)
{}


Foam::arterySteadyScalingLawInletFvPatchVectorField::
arterySteadyScalingLawInletFvPatchVectorField
(
    const arterySteadyScalingLawInletFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    area_(pivpvf.area_),
    areaAssumedDiameter_(pivpvf.areaAssumedDiameter_),
    magnitudeVelocity_(pivpvf.magnitudeVelocity_),
    values_(pivpvf.values_)
{}


Foam::arterySteadyScalingLawInletFvPatchVectorField::
arterySteadyScalingLawInletFvPatchVectorField
(
    const arterySteadyScalingLawInletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    area_(pivpvf.area_),
    areaAssumedDiameter_(pivpvf.areaAssumedDiameter_),
    magnitudeVelocity_(pivpvf.magnitudeVelocity_),
    values_(pivpvf.values_)
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
    fvPatchVectorField::write(os);
    os.writeKeyword("area") << area_ << token::END_STATEMENT << nl;
    os.writeKeyword("areaAssumedDiameter") << areaAssumedDiameter_ << token::END_STATEMENT << nl;
    os.writeKeyword("magnitudeVelocity") << magnitudeVelocity_ << token::END_STATEMENT << nl;
    values_.writeEntry("values", os);

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
