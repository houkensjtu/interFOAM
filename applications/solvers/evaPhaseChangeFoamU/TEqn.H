{
	// mass transfer
	const volScalarField& mDotT = twoPhaseProperties->mDotT();

	//mixture properties 
	volScalarField K = alpha1*K1 + (scalar(1)-alpha1)*K2;
	volScalarField rhoC = alpha1*rho1*C1 + (scalar(1) - alpha1)*rho2*C2;
	surfaceScalarField rhoPhiC("rhoPhiC", phi*fvc::interpolate(rhoC));

	//energy equation
	fvScalarMatrix TEqn
	(
		fvm::ddt(rhoC , T)
		+ fvm::div(rhoPhiC, T)
		- fvm::laplacian(K , T)
// source term for boiling 
        +L*(fvm::Sp(mDotT,T) - mDotT*TSat)

	);

	TEqn.relax();
	TEqn.solve();

	Info << "TAve = "
	<< T.weightedAverage(mesh.V()).value()
	<< "  Min(T) = " << min(T).value()
	<< "  Max(T) = " << max(T).value()
	<< endl;

}

