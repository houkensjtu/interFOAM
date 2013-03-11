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

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(phaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixture::phaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    twoPhaseMixtureI(U, phi, alpha1Name),//add
    phaseChangeTwoPhaseMixtureCoeffs_(subDict(type + "Coeffs")),
    pSat_(lookup("pSat")),
    TSat_(lookup("TSat")),//next we should consider T base on p saturation ?
    L_(lookup("lHeat"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::volScalarField Foam::phaseChangeTwoPhaseMixture::vDotAlphal() const
{

	volScalarField alphalCoeff = 1.0/rho1() - alpha1_*(1.0/rho1() - 1.0/rho2());
    volScalarField mDotAlphal = this->mDotAlphal();

    return volScalarField
    (
    		alphalCoeff*mDotAlphal
    );
}
Foam::volScalarField Foam::phaseChangeTwoPhaseMixture::vDotP() const
{
	  // dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
	//   volScalarField pCoeff = alpha1_/rho1() + (1-alpha1_)/rho2();
	   volScalarField pCoeff = 1.0/rho1() - alpha1_*(1.0/rho1() - 1.0/rho2());
	   volScalarField mDotP = this->mDotP();

	    return volScalarField
	    (
	    		pCoeff*mDotP
	    );

}

Foam::volScalarField Foam::phaseChangeTwoPhaseMixture::vDotT() const
{
	   volScalarField TCoeff = 1.0/rho2() - alpha1_*(1.0/rho2() - 1.0/rho1());
	   volScalarField mDotT = this->mDotT();

	    return volScalarField
	    (
	    		TCoeff*mDotT
	    );

}



bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (twoPhaseMixtureI::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        lookup("pSat") >> pSat_;
        lookup("TSat") >> TSat_;
        lookup("lHeat") >> L_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
