    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
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

    //-Add temperature
    Info<< "Reading field T\n"<<endl;
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

    Info<< "Reading field gradT\n"<<endl;
    volVectorField gradT
    (
      IOobject
      (
       "gradT",
       runTime.timeName(),
       mesh,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
      ),
      mesh
    );

    //-Add temperature

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi);

    volScalarField& alpha1(twoPhaseProperties.alpha1());

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

    //-Add temperature
    const dimensionedScalar& cp1 = twoPhaseProperties.cp1();
    const dimensionedScalar& cp2 = twoPhaseProperties.cp2();    
    //-Add temperature

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
            IOobject::AUTO_WRITE
        ),
        rho2*phi
    );

    //-Add temperature
    Info<<"Reading/calculating rho*cp\n"<<endl;
    volScalarField rhoCp
    (
    IOobject
    (
     "rho*Cp",
     runTime.timeName(),
     mesh,
     IOobject::NO_READ,
     IOobject::AUTO_WRITE
    ),
    alpha1*rho1*cp1+(scalar(1)-alpha1)*rho2*cp2,
    alpha1.boundaryField().types()
    );
    rhoCp.oldTime();

    Info<<"Reading/calculating rho*phi*cp\n"<<endl;
    surfaceScalarField rhoPhiCpf
    (
    IOobject
    (
     "rho*phi*cpf",
     runTime.timeName(),
     mesh,
     IOobject::NO_READ,
     IOobject::NO_WRITE
    ),
    rhoPhi*cp1
    );
    //-Add temperature

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    #include "readGravitationalAcceleration.H"

    /*
    dimensionedVector g0(g);

    // Read the data file and initialise the interpolation table
    interpolationTable<vector> timeSeriesAcceleration
    (
        runTime.path()/runTime.caseConstant()/"acceleration.dat"
    );
    */

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
