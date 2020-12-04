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

#include "pulsatilePipeInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pulsatilePipeInletFvPatchVectorField::
pulsatilePipeInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowRate_(0.),
    cardiacCycle_(1.0),
    areaTotal_(sum(patch().magSf()))
{}


Foam::pulsatilePipeInletFvPatchVectorField::
pulsatilePipeInletFvPatchVectorField
(
    const pulsatilePipeInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowRate_(0.),
    cardiacCycle_(1.0),
    areaTotal_(sum(patch().magSf()))
{
}


Foam::pulsatilePipeInletFvPatchVectorField::
pulsatilePipeInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowRate_(0.),
    cardiacCycle_(1.0),
    areaTotal_(sum(patch().magSf()))
{
    // Note: The following is the original code in this section 
    fvPatchVectorField::operator=(patch().nf());
    dict.lookup("flowRate") >> flowRate_;
    dict.lookup("cardiacCycle") >> cardiacCycle_;
    updateCoeffs();
}


Foam::pulsatilePipeInletFvPatchVectorField::
pulsatilePipeInletFvPatchVectorField
(
    const pulsatilePipeInletFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    flowRate_(pivpvf.flowRate_),
    cardiacCycle_(pivpvf.cardiacCycle_),
    areaTotal_(pivpvf.areaTotal_)
{}


Foam::pulsatilePipeInletFvPatchVectorField::
pulsatilePipeInletFvPatchVectorField
(
    const pulsatilePipeInletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    flowRate_(pivpvf.flowRate_),
    cardiacCycle_(pivpvf.cardiacCycle_),
    areaTotal_(pivpvf.areaTotal_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// extract the normals and the total area of the inlet patch

void Foam::pulsatilePipeInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::pulsatilePipeInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

// NOTE: here is the key method which implements the actual maths for calculating the
// inlet profiles
void Foam::pulsatilePipeInletFvPatchVectorField::updateCoeffs()
{
    scalar t=this->db().time().value();
    vectorField & field = *this;
    const vectorField normals = patch().nf();
    scalar coeff(
        (
            flowRate_-

            ((
                2.83586*std::cos((1*2*PI/cardiacCycle_)*t)+
                1.13367*std::sin((1*2*PI/cardiacCycle_)*t)-
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
            )/6e7)
        )/areaTotal_
    );

    if (updated())
    {
        return;
    }
    field = normals * coeff;
    // call the base class to ensure all the appropriate variables get updated
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::pulsatilePipeInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    // flowRate_.writeEntry("refValue", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pulsatilePipeInletFvPatchVectorField
    );
}

// ************************************************************************* //
