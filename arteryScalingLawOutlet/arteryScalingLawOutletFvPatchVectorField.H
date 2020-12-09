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
    Foam::arteryScalingLawOutletFvPatchVectorField

Description
    Describes a surface normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors

SourceFiles
    arteryScalingLawOutletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef arteryScalingLawOutletFvPatchVectorField_H
#define arteryScalingLawOutletFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "surfaceFields.H"
#include "cmath"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class surfaceNormalFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class arteryScalingLawOutletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    word inletPatchName_;
    word oppositeOutletPatchName_;
    scalar magnitudeVelocity_;
    vectorField values_;


public:

    //- Runtime type information
    TypeName("arteryScalingLawOutlet");


    // Constructors

        //- Construct from patch and internal field
        arteryScalingLawOutletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        arteryScalingLawOutletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  arteryScalingLawOutletFvPatchVectorField
        //  onto a new patch
        arteryScalingLawOutletFvPatchVectorField
        (
            const arteryScalingLawOutletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        arteryScalingLawOutletFvPatchVectorField
        (
            const arteryScalingLawOutletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new arteryScalingLawOutletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        arteryScalingLawOutletFvPatchVectorField
        (
            const arteryScalingLawOutletFvPatchVectorField&,
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
                new arteryScalingLawOutletFvPatchVectorField
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
