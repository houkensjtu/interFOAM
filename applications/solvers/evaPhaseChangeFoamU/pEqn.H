{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    adjustPhi(phiU, U, p);

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)*mesh.magSf()
          + fvc::interpolate(rho)*(g & mesh.Sf())
        )*rUAf;

//source term should be add in piso algorithm due to mass transfer locally 
    const volScalarField& vDotP = twoPhaseProperties-> vDotP();//add


    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUAf, p) == fvc::div(phi)
            +vDotP//add
        );


        pEqn.setReference(pRefCell, pRefValue);

        if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
        {
            pEqn.solve(mesh.solver(p.name() + "Final"));
        }
        else
        {
            pEqn.solve(mesh.solver(p.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}
