/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    vorticityTransport



Description
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicBlastFvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "wedgeFvPatch.H"
#include "compressibleSystem.H"
#include "timeIntegrator.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "addFunctionObjectOptions.H"

    #include "setRootCase.H"

    #include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"


    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        if (mesh.readUpdate() != polyMesh::UNCHANGED)
        {
            
        }
             #include "createFields.H"
              
            volVectorField vorticity
            (
                IOobject
                (
                    "vorticity",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::curl(U)
            );
            scalar fb=0.1;
            volScalarField  vorticity_dialation = mag(vorticity*fvc::div(phiByRho));
            volScalarField  vorticity_baroclinic = 1.0/(rho*rho)*mag(fvc::grad(rho)^fvc::grad(p));
            volScalarField  vorticity_stretching = mag(vorticity&fvc::grad(U)); 
            volVectorField vorticity_viscous 
            (
                IOobject
                (
                    "vorticity_viscous",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                     muEff/rho*fvc::laplacian(vorticity) // term 1  
                    - muEff/(rho*rho)*(fvc::grad(rho)^fvc::laplacian((U)))  //term 2
                    - (1.0/3.0+fb)*muEff/(rho*rho)*(fvc::grad(rho)^fvc::grad(fvc::div(phiByRho))) //term 3                        
                    +fvc::curl(muEff/rho*fvc::laplacian(U))                      //fourth term
                    + fvc::curl(2.0/rho*fvc::grad(muEff)&fvc::grad(U))      //fifth term
                    + fvc::curl(1.0/rho*fvc::grad(muEff)^vorticity)         //sixth term
                    + fvc::curl(muEff/rho*fvc::grad(fvc::div(phiByRho)))          //seventh term
                    + fvc::curl(1/rho*(fb-2.0/3.0)*fvc::grad(muEff*fvc::div(phiByRho)))
                );   //eight term
            vorticity.write();
            vorticity_viscous.write();
        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
