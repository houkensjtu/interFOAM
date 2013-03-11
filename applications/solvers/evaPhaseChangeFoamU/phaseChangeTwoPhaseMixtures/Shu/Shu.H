/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU Generac License as published by the
    Free Software Foundation; either  2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the ho it will be useful, but WITHOUT
    ANY WARRANTY; without even the imarranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.he GNU General Public License
    for more details.

    You should have received a copy oNU General Public License
    along with OpenFOAM; if not, write Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floorn, MA 02110-1301 USA

Class
    Foam::phaseChangeTwoPhaseMixtures::Shu

Description
    Shu evaporation model.

    Reference:
   // mass flux
   mflux = k grad(T)/L & n
  // mass per volume
   mdot = mflux * Ai/Vi
SourceFiles
    Shu.C

\*--------------------------------------------------------------------*/

#ifndef Shu_H
#define Shu_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Shu
\*--------------------------------------------------------------------*/

class Shu
:
    public phaseChangeTwoPhaseMixture
{
    // Private data


public:

    //- Runtime type information
    TypeName("Shu");


    // Constructors

        //- construct from components
        Shu
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Destructor

        virtual ~Shu()
        {}


    // Member Functions

        //- Return the mass vaporisation rates
        virtual volScalarField mDotAlphal() const;

        //- Correct the Shu phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
