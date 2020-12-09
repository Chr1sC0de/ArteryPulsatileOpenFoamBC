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

Class
    Foam::arteryPulsatileScalingLawInletFvPatchVectorField

Description
    Describes a surface normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors

SourceFiles
    arteryPulsatileScalingLawInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef arteryPulsatileScalingLawInletFvPatchVectorField_H
#define arteryPulsatileScalingLawInletFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "cmath"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class surfaceNormalFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class arteryPulsatileScalingLawInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
        // the area is in metres calculated by the constructors
        scalar area_;
        // area assumed diameter
        scalar areaAssumedDiameter_;
        // the flow rate is in m^3/s
        scalar meanVelocity_;
        // cardiac cycle
        scalar cardiacCycle_;
        // velocity magnitude
        scalar magnitudeVelocity_;
        // switchOn
        scalar switchOn_;
        // values
        vectorField values_;



public:

    //- Runtime type information
    TypeName("arteryPulsatileScalingLawInlet");


    // Constructors

        //- Construct from patch and internal field
        arteryPulsatileScalingLawInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        arteryPulsatileScalingLawInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  arteryPulsatileScalingLawInletFvPatchVectorField
        //  onto a new patch
        arteryPulsatileScalingLawInletFvPatchVectorField
        (
            const arteryPulsatileScalingLawInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        arteryPulsatileScalingLawInletFvPatchVectorField
        (
            const arteryPulsatileScalingLawInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new arteryPulsatileScalingLawInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        arteryPulsatileScalingLawInletFvPatchVectorField
        (
            const arteryPulsatileScalingLawInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new arteryPulsatileScalingLawInletFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }



    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );
        // Evaluation functions
        virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
