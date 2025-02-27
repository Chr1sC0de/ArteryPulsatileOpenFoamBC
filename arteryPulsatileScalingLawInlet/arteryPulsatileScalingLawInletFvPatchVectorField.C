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

#include "arteryPulsatileScalingLawInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::arteryPulsatileScalingLawInletFvPatchVectorField::
arteryPulsatileScalingLawInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    area_(0.),
    areaAssumedDiameter_(0.),
    meanVelocity_(0.),
    cardiacCycle_(0.),
    magnitudeVelocity_(0.),
    switchOn_(0.)
{}

Foam::arteryPulsatileScalingLawInletFvPatchVectorField::
arteryPulsatileScalingLawInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    area_(0.),
    areaAssumedDiameter_(0.),
    meanVelocity_(0.),
    cardiacCycle_(0.),
    magnitudeVelocity_(0.),
    switchOn_(0.)
{
    // Note: The following is the original code in this section
    // fvPatchVectorField::operator=(patch().nf());
    area_                = dict.lookupOrDefault<scalar>("area", sum(patch().magSf()));
    areaAssumedDiameter_ = dict.lookupOrDefault<scalar>("areaAssumedDiameter", 2*std::sqrt(area_/PI));
    meanVelocity_        = dict.lookupOrDefault<scalar>("meanVelocity", 1.43*std::pow(areaAssumedDiameter_, 2.55)/area_);
    cardiacCycle_        = dict.lookupOrDefault<scalar>("cardiacCycle", 0.8);
    switchOn_            = dict.lookupOrDefault<scalar>("switchOn", 0.);

    if (dict.found("magnitudeVelocity")){
        dict.lookup("magnitudeVelocity") >> magnitudeVelocity_;
        fvPatchVectorField::operator=(-patch().nf() * magnitudeVelocity_);
    } else {
        updateCoeffs();
    }
}

Foam::arteryPulsatileScalingLawInletFvPatchVectorField::
arteryPulsatileScalingLawInletFvPatchVectorField
(
    const arteryPulsatileScalingLawInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    area_(ptf.area_),
    areaAssumedDiameter_(ptf.areaAssumedDiameter_),
    meanVelocity_(ptf.meanVelocity_),
    cardiacCycle_(ptf.cardiacCycle_),
    magnitudeVelocity_(ptf.magnitudeVelocity_),
    switchOn_(ptf.switchOn_)
{}


Foam::arteryPulsatileScalingLawInletFvPatchVectorField::
arteryPulsatileScalingLawInletFvPatchVectorField
(
    const arteryPulsatileScalingLawInletFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    area_(pivpvf.area_),
    areaAssumedDiameter_(pivpvf.areaAssumedDiameter_),
    meanVelocity_(pivpvf.meanVelocity_),
    cardiacCycle_(pivpvf.cardiacCycle_),
    magnitudeVelocity_(pivpvf.magnitudeVelocity_),
    switchOn_(pivpvf.switchOn_)
{}


Foam::arteryPulsatileScalingLawInletFvPatchVectorField::
arteryPulsatileScalingLawInletFvPatchVectorField
(
    const arteryPulsatileScalingLawInletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    area_(pivpvf.area_),
    areaAssumedDiameter_(pivpvf.areaAssumedDiameter_),
    meanVelocity_(pivpvf.meanVelocity_),
    cardiacCycle_(pivpvf.cardiacCycle_),
    magnitudeVelocity_(pivpvf.magnitudeVelocity_),
    switchOn_(pivpvf.switchOn_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// extract the normals and the total area of the inlet patch

void Foam::arteryPulsatileScalingLawInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::arteryPulsatileScalingLawInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

// NOTE: here is the key method which implements the actual maths for calculating the
// inlet profiles
void Foam::arteryPulsatileScalingLawInletFvPatchVectorField::updateCoeffs()
{

    if (updated())
    {
        return;
    }

    const vectorField normals = patch().nf();
    scalar t                  = this->db().time().value();

    if (t >= switchOn_){

        t -= switchOn_;

        magnitudeVelocity_=
            meanVelocity_-
            (
                2.83586*std::cos((1*2*PI/cardiacCycle_)*t)+
                0.13367*std::sin((1*2*PI/cardiacCycle_)*t)-
                1.09937*std::cos((2*2*PI/cardiacCycle_)*t)-
                0.85915*std::sin((2*2*PI/cardiacCycle_)*t)+
                0.70761*std::cos((3*2*PI/cardiacCycle_)*t)-
                1.04898*std::sin((3*2*PI/cardiacCycle_)*t)+
                0.48890*std::cos((4*2*PI/cardiacCycle_)*t)+
                0.71584*std::sin((4*2*PI/cardiacCycle_)*t)-
                0.29543*std::cos((5*2*PI/cardiacCycle_)*t)-
                0.06778*std::sin((5*2*PI/cardiacCycle_)*t)+
                0.11454*std::cos((6*2*PI/cardiacCycle_)*t)+
                0.22741*std::sin((6*2*PI/cardiacCycle_)*t)-
                0.37011*std::cos((7*2*PI/cardiacCycle_)*t)-
                0.13848*std::sin((7*2*PI/cardiacCycle_)*t)+
                0.23613*std::cos((8*2*PI/cardiacCycle_)*t)-
                0.15284*std::sin((8*2*PI/cardiacCycle_)*t)-
                0.04670*std::cos((9*2*PI/cardiacCycle_)*t)+
                0.08242*std::sin((9*2*PI/cardiacCycle_)*t)+
                0.08075*std::cos((10*2*PI/cardiacCycle_)*t)-
                0.07041*std::sin((10*2*PI/cardiacCycle_)*t)+
                0.00550*std::cos((11*2*PI/cardiacCycle_)*t)+
                0.12147*std::sin((11*2*PI/cardiacCycle_)*t)-
                0.07736*std::cos((12*2*PI/cardiacCycle_)*t)-
                0.03547*std::sin((12*2*PI/cardiacCycle_)*t)+
                0.04093*std::cos((13*2*PI/cardiacCycle_)*t)-
                0.01490*std::sin((13*2*PI/cardiacCycle_)*t)-
                0.03809*std::cos((14*2*PI/cardiacCycle_)*t)+
                0.01484*std::sin((14*2*PI/cardiacCycle_)*t)+
                0.03868*std::cos((15*2*PI/cardiacCycle_)*t)-
                0.05287*std::sin((15*2*PI/cardiacCycle_)*t)
            )/6e7/area_ ;
    } else {
        magnitudeVelocity_ = meanVelocity_;
    }
    this->operator == (-normals * magnitudeVelocity_);
    // call the base class to ensure all the appropriate variables get updated
    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::arteryPulsatileScalingLawInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("area") << area_ << token::END_STATEMENT << nl;
    os.writeKeyword("areaAssumedDiameter") << areaAssumedDiameter_ << token::END_STATEMENT << nl;
    os.writeKeyword("meanVelocity") << meanVelocity_ << token::END_STATEMENT << nl;
    os.writeKeyword("cardiacCycle") << cardiacCycle_ << token::END_STATEMENT << nl;
    os.writeKeyword("magnitudeVelocity") << magnitudeVelocity_ << token::END_STATEMENT << nl;
    os.writeKeyword("switchOn") << switchOn_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        arteryPulsatileScalingLawInletFvPatchVectorField
    );
}

// ************************************************************************* //
