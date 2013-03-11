{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();


    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            );

  //source terms in void fraction equation
   const volScalarField& vDotAlphal = twoPhaseProperties-> vDotAlphal();//add


        volScalarField Sp //add
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
         - vDotAlphal

        );

        volScalarField Su //add
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
    
         mesh,
         dimensionedScalar("Su",dimensionSet(0,0,-1,0,0,0,0),0)
        );

        MULES::implicitSolve(geometricOneField(), alpha1, phi, phiAlpha, Sp, Su, 1, 0);

//        MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);
    //    MULES::explicitSolve(oneField(), alpha1, phi, phiAlpha, Sp, Su, 1, 0);
        rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
         << endl;
}
