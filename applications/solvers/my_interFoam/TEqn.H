volScalarField kappa = twoPhaseProperties.kappa();
volScalarField Cp = alpha1*cp1 + (scalar(1) - alpha1)*cp2;

forAll(mesh.cells(),celli)
{
  if(alpha1[celli]<scalar(0.70))
    T[celli]=373;
};
fvScalarMatrix TEqn
(
 fvm::ddt(rho,T)
 +fvm::div(rhoPhi,T)
 -fvm::laplacian(kappa/Cp,T)
);

TEqn.solve();
