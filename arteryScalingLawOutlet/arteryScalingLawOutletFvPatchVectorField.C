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

#include "arteryScalingLawOutletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::arteryScalingLawOutletFvPatchVectorField::
arteryScalingLawOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletPatchName_("INLET"),
    oppositeOutletPatchName_("OUTLET_2"),
    magnitudeVelocity_(0.)
{}


Foam::arteryScalingLawOutletFvPatchVectorField::
arteryScalingLawOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletPatchName_("INLET"),
    oppositeOutletPatchName_("OUTLET_2"),
    magnitudeVelocity_(0.)
{
    // Note: No need to have an updateCoeff method as the inlet value is constant 
    // fvPatchVectorField::operator=(patch().nf()*inletVelocity_);
    // lets start by reading in the inlet and outlet patch names
    // fvPatchVectorField::operator=(patch().nf()*0.);

    inletPatchName_          = dict.lookupOrDefault<word>("inletPatchName", "INLET");
    oppositeOutletPatchName_ = dict.lookupOrDefault<word>("oppositeOutletPatchName", "OUTLET_2");
    if (dict.found("magnitudeVelocity")){
        dict.lookup("magnitudeVelocity") >> magnitudeVelocity_;
        fvPatchVectorField::operator=(patch().nf() * magnitudeVelocity_);
    } else {
        updateCoeffs();
    }
}


Foam::arteryScalingLawOutletFvPatchVectorField::
arteryScalingLawOutletFvPatchVectorField
(
    const arteryScalingLawOutletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletPatchName_(ptf.inletPatchName_),
    oppositeOutletPatchName_(ptf.oppositeOutletPatchName_),
    magnitudeVelocity_(ptf.magnitudeVelocity_)
{}


Foam::arteryScalingLawOutletFvPatchVectorField::
arteryScalingLawOutletFvPatchVectorField
(
    const arteryScalingLawOutletFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    inletPatchName_(pivpvf.inletPatchName_),
    oppositeOutletPatchName_(pivpvf.oppositeOutletPatchName_),
    magnitudeVelocity_(pivpvf.magnitudeVelocity_)
{}


Foam::arteryScalingLawOutletFvPatchVectorField::
arteryScalingLawOutletFvPatchVectorField
(
    const arteryScalingLawOutletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    inletPatchName_(pivpvf.inletPatchName_),
    oppositeOutletPatchName_(pivpvf.oppositeOutletPatchName_),
    magnitudeVelocity_(pivpvf.magnitudeVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// extract the normals and the total area of the inlet patch

void Foam::arteryScalingLawOutletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::arteryScalingLawOutletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::arteryScalingLawOutletFvPatchVectorField::updateCoeffs()
{
    if (updated()){ return; }
    // get the boundary mesh and the surface scalar fields
    const fvMesh & mesh              = patch().boundaryMesh().mesh();
    const surfaceScalarField & magSf = mesh.magSf();
    // find the labels for the inlet and oulet patchs
    const label inletLabel         = mesh.boundaryMesh().findPatchID(inletPatchName_);
    const label oppositePatchLabel = mesh.boundaryMesh().findPatchID(oppositeOutletPatchName_);
    // get the areas of each of the patches
    const scalarField & inletAreas_ = magSf.boundaryField()[inletLabel];
    const scalar inletPatchAreaTotal_          = gSum(inletAreas_);
    const scalar currentOutletPatchAreaTotal_  = gSum(patch().magSf());
    const scalar oppositeOutletPatchAreaTotal_ = gSum(magSf.boundaryField()[oppositePatchLabel]);
    // get the area assumed diameters of each of the inlet
    const scalar inletPatchAreaAssumedDiameter_          = 2*std::sqrt(inletPatchAreaTotal_/PI);
    // get the flow rate ratio of the inlet and outlet diameters 
    const scalar q_opp_over_q_current = std::pow(oppositeOutletPatchAreaTotal_/currentOutletPatchAreaTotal_, 2.27);
    // get the current time steps and set a variable for the inlet flow rate

    scalar t = this->db().time().value();
    scalar inletPatchFlowRate_;

    // now calculate the flow rate of the inlet patch, if the time is 0 we need to assume an inlet velocity
    // as the velocity fields have not been calculated
    if (t==0.)
    {
        inletPatchFlowRate_ = 1.43*std::pow(inletPatchAreaAssumedDiameter_, 2.55);
    } else {
        scalarField magInletVelocities_ = 
            mag(mesh.boundary()[inletLabel].lookupPatchField<volVectorField, vector>("U"));
        inletPatchFlowRate_ = gSum(inletAreas_*magInletVelocities_);
    }

    scalar currentOutletPatchFlowRate_ = inletPatchFlowRate_ / (q_opp_over_q_current + 1);

    magnitudeVelocity_ = 
        currentOutletPatchFlowRate_/currentOutletPatchAreaTotal_;
    //update the velocity field of the current outlet patch
    this->operator==(patch().nf() * magnitudeVelocity_);
    // call the base class to ensure all the appropriate variables get updated
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::arteryScalingLawOutletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("inletPatchName") << inletPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("oppositeOutletPatchName") << oppositeOutletPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("magnitudeVelocity") << magnitudeVelocity_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        arteryScalingLawOutletFvPatchVectorField
    );
}

// ************************************************************************* //
