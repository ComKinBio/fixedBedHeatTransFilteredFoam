rhoc = thermo.rho();

volScalarField rAUc(1.0/UcEqn.A());
surfaceScalarField rhocrAUcf("rhorAUf", fvc::interpolate(rhoc*rAUc));
volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, p));

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            fvc::flux(HbyA)
          + MRF.zeroFilter
            (
                rhocrAUcf*fvc::ddtCorr(rhoc, Uc, rhocPhic)/fvc::interpolate(rhoc)
            )
        )
    );

    MRF.makeRelative(fvc::interpolate(psi), phid);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(alphacavg, psi, p)
          + fvm::div(alphacf*phid, p)
          - fvm::laplacian(alphacf*alphacf*rhocrAUcf, p)
         ==
            fvOptions(psi, p, rhoc.name())
        );

//         pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            rhocPhic == pEqn.flux()/alphacf;
            phic = rhocPhic/fvc::interpolate(rhoc);
        }
    }
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rhoc*HbyA)
          + MRF.zeroFilter(alphacf*rhocrAUcf*fvc::ddtCorr(rhoc, Uc, rhocPhic))
        )
    );

    if (p.needReference())
    {
        adjustPhi(phiHbyA, Uc, p);
    }

    MRF.makeRelative(fvc::interpolate(rhoc), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rhoc, Uc, phiHbyA, rhocrAUcf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(alphacavg, psi, p)
          + fvc::div(alphacf*phiHbyA)
          - fvm::laplacian(alphacf*alphacf*rhocrAUcf, p)
         ==
            fvOptions(psi, p, rhoc.name())
        );

//         pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            rhocPhic = phiHbyA + pEqn.flux()/alphacf;
            phic = rhocPhic/fvc::interpolate(rhoc);
        }
    }
}


alphaRhoPhic = alphacf*rhocPhic;

#include "rhocEqn.H"
#include "compressibleContinuityErrs.H"

Uc = HbyA - rAUc*alphacavg*fvc::grad(p);
Uc.correctBoundaryConditions();
fvOptions.correct(Uc);

Kc = 0.5*magSqr(Uc);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
