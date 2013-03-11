volScalarField kappa = twoPhaseProperties.kappa();
volScalarField Cp = alpha1*cp1 + (scalar(1) - alpha1)*cp2;
surfaceScalarField kappaf = twoPhaseProperties.kappaf();

forAll(mesh.cells(),celli)
{
  if(alpha1[celli]<scalar(0.30))
    T[celli]=373;
};

fvScalarMatrix TEqn
(
 fvm::ddt(rhoCp,T)
 +fvm::div(rhoPhiCpf,T)
 -fvm::laplacian(kappaf,T)
);

TEqn.solve();
gradT = fvc::grad(T);
