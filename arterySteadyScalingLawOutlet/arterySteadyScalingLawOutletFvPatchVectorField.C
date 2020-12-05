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

#include "arterySteadyScalingLawOutletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::arterySteadyScalingLawOutletFvPatchVectorField::
arterySteadyScalingLawOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletPatchName_("INLET"),
    oppositeOutletPatchName_("OUTLET_2")
{}


Foam::arterySteadyScalingLawOutletFvPatchVectorField::
arterySteadyScalingLawOutletFvPatchVectorField
(
    const arterySteadyScalingLawOutletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletPatchName_("INLET"),
    oppositeOutletPatchName_("OUTLET_2")
{}


Foam::arterySteadyScalingLawOutletFvPatchVectorField::
arterySteadyScalingLawOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletPatchName_("INLET"),
    oppositeOutletPatchName_("OUTLET")
{
    // Note: No need to have an updateCoeff method as the inlet value is constant 
    // fvPatchVectorField::operator=(patch().nf()*inletVelocity_);
    // lets start by reading in the inlet and outlet patch names
    inletPatchName_          = dict.lookupOrDefault<word>("inletPatchName", "INLET");
    oppositeOutletPatchName_ = dict.lookupOrDefault<word>(
        "oppositeOutletPatchName", "OUTLET");

    const fvMesh & mesh              = patch().boundaryMesh().mesh();
    const surfaceScalarField & magSf = mesh.magSf();

    const label inletLabel         = mesh.boundaryMesh().findPatchID(inletPatchName_);
    const label oppositePatchLabel = mesh.boundaryMesh().findPatchID(
        oppositeOutletPatchName_);
    // get the area of each of the patches
    const scalar currentOutletPatchAreaTotal_  = gSum(patch().magSf());
    const scalar inletPatchAreaTotal_          = gSum(magSf.boundaryField()[inletLabel]);
    const scalar oppositeOutletPatchAreaTotal_ = gSum(magSf.boundaryField()[oppositePatchLabel]);
    // get the area assumed dimaters of each patch
    const scalar currentOutletPatchAreaAssumedDiameter_  = 2*std::sqrt(
        currentOutletPatchAreaTotal_/PI);
    const scalar inletPatchAreaAssumedDiameter_          = 2*std::sqrt(inletPatchAreaTotal_/PI);
    const scalar oppositeOutletPatchAreaAssumedDiameter_ = 2*std::sqrt(
        oppositeOutletPatchAreaTotal_/PI);
    // get the flow rate of the inlet patch
    const scalar inletPatchFlowRate_ = 1.43*std::pow(inletPatchAreaAssumedDiameter_, 2.55);
    // get the outlflow ratio
    const scalar q_opp_over_q_current = std::pow(
        oppositeOutletPatchAreaTotal_/currentOutletPatchAreaTotal_, 2.27);
    // get the flow rates for the opposite oulet and the current outlet
    const scalar oppositeOutletPatchFlowRate_ = 
        q_opp_over_q_current * inletPatchFlowRate_ / (q_opp_over_q_current + 1);
    const scalar currentOutletPatchFlowRate_ = 
        inletPatchFlowRate_ - oppositeOutletPatchFlowRate_;
    // get the patch velocity 
    const scalar currenOutletPatchVelocity_ = 
        currentOutletPatchFlowRate_/currentOutletPatchAreaTotal_;

    fvPatchVectorField::operator=(patch().nf()*currenOutletPatchVelocity_);

}


Foam::arterySteadyScalingLawOutletFvPatchVectorField::
arterySteadyScalingLawOutletFvPatchVectorField
(
    const arterySteadyScalingLawOutletFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    inletPatchName_(pivpvf.inletPatchName_),
    oppositeOutletPatchName_(pivpvf.oppositeOutletPatchName_)
{}


Foam::arterySteadyScalingLawOutletFvPatchVectorField::
arterySteadyScalingLawOutletFvPatchVectorField
(
    const arterySteadyScalingLawOutletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    inletPatchName_(pivpvf.inletPatchName_),
    oppositeOutletPatchName_(pivpvf.oppositeOutletPatchName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// extract the normals and the total area of the inlet patch

void Foam::arterySteadyScalingLawOutletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::arterySteadyScalingLawOutletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void Foam::arterySteadyScalingLawOutletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("inletPatchName") << inletPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("oppositeOutletPatchName") << oppositeOutletPatchName_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        arterySteadyScalingLawOutletFvPatchVectorField
    );
}

// ************************************************************************* //
