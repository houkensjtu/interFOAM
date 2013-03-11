// mass transfer
const volScalarField& mDotT = twoPhaseProperties->mDotT();

volScalarField kappa = twoPhaseProperties->kappa(); 
volScalarField Cp = alpha1*cp1 + (scalar(1) - alpha1)*cp2; 
surfaceScalarField kappaf = twoPhaseProperties->kappaf();

//
//forAll(mesh.cells(),celli)
//{
//  if(alpha1[celli] > scalar(0.50))
//    T[celli]=373.15;
//};

fvScalarMatrix TEqn
(
 fvm::ddt(rhoCp,T)
 +fvm::div(rhoPhiCpf,T)
 -fvm::laplacian(kappaf,T)

// source term for boiling 
 +L*(fvm::Sp(mDotT,T) - mDotT*TSat)

);

TEqn.solve();
gradT = fvc::grad(T);
