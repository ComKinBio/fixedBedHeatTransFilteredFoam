/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    coalChemistryFoam

Group
    grpLagrangianSolvers

Description
    Transient solver for compressible, turbulent flow, with coal and limestone
    particle clouds, an energy source, and combustion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "coalChemistryTurbulenceModel.H"
#include "basicThermoCloud.H"
#include "coalCloud.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "fvOptions.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "simpleControl.H" //added for DBM
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "coarserGrid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

//     #include "addCheckCaseOptions.H"
//     #include "setRootCase.H"
    #include "setRootCaseLists.H" //Tian
    #include "createTime.H"
    #include "createMesh.H"
//     #include "createControl.H" //acture this file is to creat control object according to the head file, here we include two head file, it cause error using this file.
    pimpleControl pimple(mesh);
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    #include "diffusion.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        
        
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        
        coalParcels.evolve();
//         limestoneParcels.evolve();

//         alphac = max(1.0 - coalParcels.theta() /*- limestoneParcels.theta()*/, alphacMin);
//         alphac.correctBoundaryConditions();
//         
//         if (useDiffusionMethod)
//         {
//             #include "alphaDiffusionFunc.H"
//             alphac = max(alphac, alphacMin);
//             alphac.correctBoundaryConditions();
//         }
//         
//         if (useCoarseGridSource)
//         {
//             alphac = max(alphacavg, alphacMin);
//             alphac.correctBoundaryConditions();
//         }
//         
//         alphacf = fvc::interpolate(alphac);
//         //phic = linearInterpolate(Uc) & mesh.Sf();
//         //rhocPhic = fvc::interpolate(rhoc)*phic;
//         alphaRhoPhic = alphacf*rhocPhic;
//         
//         
// 
//         volScalarField::Internal averagedSrho = coalParcels.Srho();
//         fvScalarMatrix averagedSrhoRho = coalParcels.Srho(rhoc);
//         fvVectorMatrix averagedSU = coalParcels.SU(Uc);
//         
//         
//         #include "rhocEqn.H"
// 
//         // --- Pressure-velocity PIMPLE corrector loop
//         while (pimple.loop())
//         {
//             #include "UcEqn.H"
//             
//             // --- Pressure corrector loop
//             while (pimple.correct())
//             {
//                 #include "pEqn.H"
//             }
// 
//             if (pimple.turbCorr())
//             {
//                 turbulence->correct();
//             }
//         }
// 
//         rhoc = thermo.rho();

        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
