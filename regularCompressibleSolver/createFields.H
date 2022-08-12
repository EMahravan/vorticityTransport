            autoPtr<fluidThermo> pThermo
            (
                fluidThermo::New(mesh)
            );
            fluidThermo& thermo = pThermo();
            thermo.validate(args.executable(), "h", "e");

            volScalarField& p = thermo.p();

            volScalarField rho
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                thermo.rho()
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
            #include "compressibleCreatePhi.H"
            Info<< "Creating turbulence model\n" << endl;
            autoPtr<compressible::momentumTransportModel> turbulence
            (
                compressible::momentumTransportModel::New
                (
                    rho,
                    U,
                    phi,
                    thermo
                )
            );

            Info<< "Creating thermophysical transport model\n" << endl;
            autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
            (
                fluidThermophysicalTransportModel::New(turbulence(), thermo)
            );

            dimensionedScalar initialMass = fvc::domainIntegrate(rho);
            volScalarField muEff("muEff", turbulence->muEff());
            
surfaceScalarField phiByRho
(
    IOobject
    (
        "phiByRho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);
