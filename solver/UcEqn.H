    MRF.correctBoundaryVelocity(Uc);

    fvVectorMatrix UcEqn(Uc, rhoc.dimensions()*Uc.dimensions()*dimVol/dimTime);

    UcEqn = 
    (
        fvm::ddt(alphacavg, rhoc, Uc) + fvm::div(alphaRhoPhic, Uc)
      + MRF.DDt(alphacavg*rhoc, Uc)
      + turbulence->divDevRhoReff(Uc)
     ==
        alphacavg*rhoc()*g
      + averagedSU
      + fvOptions(rhoc, Uc)
    );
// 
    UcEqn.relax();

    fvOptions.constrain(UcEqn);

    if (pimple.momentumPredictor())
    {
        solve(UcEqn == -alphacavg*fvc::grad(p));

        fvOptions.correct(Uc);
        Kc = 0.5*magSqr(Uc);
    }
