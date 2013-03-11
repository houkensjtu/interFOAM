    surfaceScalarField muEff
    (
        "muEff",
      //  twoPhaseProperties.muf()
        twoPhaseProperties->muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)//add
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rho)*(g & mesh.Sf())
              + (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - fvc::snGrad(p)
                ) * mesh.magSf()
            )
        );
    }
