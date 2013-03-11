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

Class
    Foam::twoPhaseThermalMixture

Description
    Extends the twoPhaseMixture class (a two-phase incompressible 
    transportModel) by adding data members and functions for fluid 
    properties relevant to energy transport. 
    
    Written by: 
	Emmanuel Doro
    Georgia Institute of Technology
    Atlanta, USA
    January 2010.

SourceFiles
    twoPhaseThermalMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseThermalMixture_H
#define twoPhaseThermalMixture_H

#include "incompressible/incompressibleTwoPhaseMixture/twoPhaseMixture.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseThermalMixture Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseThermalMixture
:
    public twoPhaseMixture
{
protected:

    // Protected data

        // specific heats
        dimensionedScalar cp1_;
        dimensionedScalar cp2_;
        
        // Prandtl numbers
        dimensionedScalar Pr1_;
        dimensionedScalar Pr2_;

        // thermal conductivity
        volScalarField kappa_;


    // Private Member Functions

        // calculates thermal conductivity
        void calcKappa();


public:

    // Constructors

        //- Construct from components
        twoPhaseThermalMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Destructor

        ~twoPhaseThermalMixture()
        {}


    // Member Functions

        //- Return const-access to phase1 specific heat
		const dimensionedScalar& cp1() const
		{
    		return cp1_;
		}
		
		//- Return const-access to phase2 specific heat
		const dimensionedScalar& cp2() const
		{
    		return cp2_;
		};
		
		//- Return const-access to phase1 Prandtl number
		const dimensionedScalar& Pr1() const
		{
    		return Pr1_;
		}
		
		//- Return const-access to phase2 Prandtl number
		const dimensionedScalar& Pr2() const
		{
    		return Pr2_;
		};

        //- Return thermal conductivity
        tmp<volScalarField> kappa() const
        {
        	return kappa_;
        }
        
        //- Return the face-interpolated thermal conductivity
        tmp<surfaceScalarField> kappaf() const;
        
        //- Return thermal diffusivity
        tmp<volScalarField> thermalDiffusivity() const;
        
        //- Return face-interpolated thermal diffusivity
        tmp<surfaceScalarField> thermalDiffusivityFace() const;

        //- Correct transport properties
        virtual void correct()
        {
            twoPhaseMixture::correct();
            calcKappa();        
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
