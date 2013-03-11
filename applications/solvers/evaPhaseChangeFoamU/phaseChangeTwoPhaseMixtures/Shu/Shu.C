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
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Shu.H"
#include "fvc.H"//add
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Shu, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeTwoPhaseMixture,
        Shu,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Shu::Shu
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi, alpha1Name)

{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam :: volScalarField Foam::phaseChangeTwoPhaseMixtures::Shu::mDotAlphal() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");

    volScalarField	Tg =T*(1-alpha1_)+TSat_*alpha1_;
    volVectorField gradAlpha = fvc::grad(alpha1_);

    volScalarField K=K1_*alpha1_+K2_*(1-alpha1_);


    return volScalarField
    (
    		scalar(-1.0)*(K/L_)*(fvc::grad(Tg)& gradAlpha)

    );

}

void Foam::phaseChangeTwoPhaseMixtures::Shu::correct()
{}


bool Foam::phaseChangeTwoPhaseMixtures::Shu::read()
{

    if (phaseChangeTwoPhaseMixture::read())
    {
        //phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type());
    	//+ "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
