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

#include "Hardt.H"
#include "fvc.H"//add
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Hardt, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeTwoPhaseMixture,
        Hardt,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Hardt::Hardt
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi, alpha1Name),

    R_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("R")),
    Cv_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv")),
    Tc_("Tc",dimensionSet(0,0,0,1,0,0,0),scalar(0.0))// Temperature criteria
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam :: dimensionedScalar Foam::phaseChangeTwoPhaseMixtures::Hardt::Rint() const
{
	return (
			(2.0-Cv_)*pow(2.0*3.1415*R_*pow(TSat_,3.0),0.5)/
			(2.0*Cv_*pow(L_,2)*rho2_)
	);
}
Foam ::volScalarField Foam::phaseChangeTwoPhaseMixtures::Hardt::intV(
		) const
{
	const scalarField vol = alpha1_.mesh().V();

	volScalarField limitedAlpha1 = alpha1_;
	
	volVectorField gradAlpha = fvc::grad(limitedAlpha1);
	volScalarField alphaPrim = mag(gradAlpha);

	dimensionedScalar n = 1.0;
	//check possible get inf
	if (sum(alphaPrim*vol*limitedAlpha1)!= 0.0)
	n=sum(alphaPrim*vol)/sum(alphaPrim*vol*limitedAlpha1);

	return(
			alphaPrim

			);
}

Foam :: volScalarField Foam::phaseChangeTwoPhaseMixtures::Hardt::mDotAlphal() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    volScalarField intV = this->intV();

return volScalarField
    (
    	max((T-TSat_),Tc_) *intV/(Rint()*L_)
    );

}

Foam :: volScalarField Foam::phaseChangeTwoPhaseMixtures::Hardt::mDotP() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    volScalarField intV = this->intV();
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    return volScalarField
    (
    		 max((T-TSat_),Tc_)
    		*intV/(Rint()*L_)
    		*limitedAlpha1
    );
}

Foam :: volScalarField Foam::phaseChangeTwoPhaseMixtures::Hardt::mDotT() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    volScalarField intV = this->intV();
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    return volScalarField
    (
    		pos((T-TSat_))
    		*intV/(Rint()*L_)
    		*limitedAlpha1
    );
}



void Foam::phaseChangeTwoPhaseMixtures::Hardt::correct()
{}


bool Foam::phaseChangeTwoPhaseMixtures::Hardt::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("R") >> R_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv") >> Cv_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
