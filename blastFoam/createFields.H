
// Read gravitational acceleration
Info<< "\nReading g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedVector(dimAcceleration, Zero)
);

autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
autoPtr<compressibleSystem> fluid
(
    compressibleSystem::New(mesh)
);
integrator->addSystem(fluid());

const volScalarField& p = fluid->p();
const volScalarField& rho = fluid->rho();
const volVectorField& U = fluid->U();

fluid->update();

volScalarField muEff("muEff", fluid->turbulence().muEff());
            
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
