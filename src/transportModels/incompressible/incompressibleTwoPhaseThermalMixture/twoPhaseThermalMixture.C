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
	
	
	Written by: 
	Emmanuel Doro
    Georgia Institute of Technology
    Atlanta, USA
    January 2010.
\*---------------------------------------------------------------------------*/

#include "twoPhaseThermalMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate thermal conductivity
void twoPhaseThermalMixture::calcKappa()
{
	volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    kappa_ = limitedAlpha1*rho1_*cp1_*(1.0/Pr1_)*nuModel1_->nu() 
    		 + (scalar(1) - limitedAlpha1)*rho2_*cp2_*(1.0/Pr2_)*nuModel2_->nu();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhaseThermalMixture::twoPhaseThermalMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    twoPhaseMixture(U, phi, alpha1Name),
        
    cp1_(nuModel1_->viscosityProperties().lookup("cp")),
	cp2_(nuModel2_->viscosityProperties().lookup("cp")),
	Pr1_(nuModel1_->viscosityProperties().lookup("Pr")),
	Pr2_(nuModel2_->viscosityProperties().lookup("Pr")),

    kappa_
    (
        IOobject
        (
            "kappa",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("kappa", dimensionSet(1, 1, -3, -1, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    calcKappa();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> twoPhaseThermalMixture::kappaf() const
{
    surfaceScalarField alpha1f =
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1));
        
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "kappaf",
            alpha1f*rho1_*cp1_*(1.0/Pr1_)*fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - alpha1f)*rho2_*cp2_*(1.0/Pr2_)*fvc::interpolate(nuModel2_->nu())
        )
    );
}


tmp<volScalarField> twoPhaseThermalMixture::thermalDiffusivity() const
{
	volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
	
	return tmp<volScalarField>
	(
    	new volScalarField
    	(
       		"thermalAlpha",
        	kappa()/(limitedAlpha1*rho1_*cp1_ +	(scalar(1) - limitedAlpha1)*rho2_*cp2_)
        )
	);
}


tmp<surfaceScalarField> twoPhaseThermalMixture::thermalDiffusivityFace() const
{
	surfaceScalarField alpha1f =
		min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1));
	
	return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "thermalAphaf",
            kappaf()/(alpha1f*rho1_*cp1_ +	(scalar(1) - alpha1f)*rho2_*cp2_)
        )
    );

}


bool twoPhaseThermalMixture::read()
{
    if (twoPhaseMixture::read())
    {
        nuModel1_->viscosityProperties().lookup("cp") >> cp1_;
        nuModel2_->viscosityProperties().lookup("cp") >> cp2_;
        nuModel1_->viscosityProperties().lookup("Pr") >> Pr1_;
        nuModel2_->viscosityProperties().lookup("Pr") >> Pr2_;
        
        return true;
     }
     else
     {
     	return false;
     }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
