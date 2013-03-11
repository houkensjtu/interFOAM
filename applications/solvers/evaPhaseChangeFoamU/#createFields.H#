    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;//add T
    volScalarField T
    (
    		IOobject
    		(
    			"T",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
    		),
    		mesh
    );


    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
    		IOobject
    		(
    			"alpha1",
    			runTime.timeName(),
    			mesh,
    			IOobject::MUST_READ,
    			IOobject::AUTO_WRITE
			),
			mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
   //twoPhaseMixtureI twoPhaseProperties(U, phi); //add twoPhaseMixtureI
    autoPtr<phaseChangeTwoPhaseMixture> twoPhaseProperties =
    		phaseChangeTwoPhaseMixture::New(U, phi);

    const dimensionedScalar& rho1 = twoPhaseProperties->rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties->rho2();

    const dimensionedScalar& K1 = twoPhaseProperties->K1();
    const dimensionedScalar& K2 = twoPhaseProperties->K2();

    const dimensionedScalar& C1 = twoPhaseProperties->C1();
    const dimensionedScalar& C2 = twoPhaseProperties->C2();

    const dimensionedScalar& L = twoPhaseProperties->lHeat();
    const dimensionedScalar& TSat = twoPhaseProperties->TSat();

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

    //A more suitable distribution of source term
    //can be computed by solving an inhomogeneous Helmholtz equation for an additional scalarField
    //helmHoltz sDis(alpha1, twoPhaseProperties); //add


     // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties())//why?
    );

    //*********would be earased ************************
    /*
        volScalarField mDotCorr //add
        (
            IOobject
            (
                "mDotCorr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar ("mDotCorr",dimensionSet(1,-3,-1,0,0,0,0),0.0),
            zeroGradientFvPatchScalarField::typeName
        );
    */
        volScalarField mDot //add
            (
                IOobject
                (
                    "mDot",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar ("mDot",dimensionSet(1,-3,-1,0,0,0,0),0.0)
               // zeroGradientFvPatchScalarField::typeName
            );

        volScalarField vDot //add
            (
                IOobject
                (
                    "vDot",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar ("vDot",dimensionSet(0,0,-1,0,0,0,0),0.0)
               // zeroGradientFvPatchScalarField::typeName
            );

        volScalarField gradAlpha //add
        (
            IOobject
            (
                "gradAlpha",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mag(fvc::grad(alpha1))
        );
