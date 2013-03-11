surfaceScalarField kappa = twoPhaseProperties.kappa();
surfaceScalarField Cp = alpha1*cp1 + (scalar(1) - alpha1)*cp2;

fvScalarMatrix TEqn
(
 fvm::ddt(rho,T)
 +fvm::div(rhoPhi,T)
 -fvm::laplacian((1.0/Cp),T)
);

TEqn.solve();
