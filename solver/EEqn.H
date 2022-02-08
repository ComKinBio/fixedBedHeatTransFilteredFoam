{
    volScalarField& he = thermo.he();

    const volScalarField alphaEff(turbulence->alphaEff());
    
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            alphaRhoPhic,
            mesh.divScheme("div(alphaRhoPhi,he)")
        )
    );

    fvScalarMatrix EEqn
    (
        fvm::ddt(alphacavg, rhoc, he) + mvConvection->fvmDiv(alphaRhoPhic, he)
      + fvc::ddt(alphacavg, rhoc, Kc) + fvc::div(alphaRhoPhic, Kc)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(alphaRhoPhic/fvc::interpolate(rhoc), Uc),
                p,
                "div(phiv,p)"
            )
          : -dpdt*alphacavg
        )
      - fvm::laplacian                
        (
            fvc::interpolate(alphacavg)
           *fvc::interpolate(alphaEff),
            he,
            "laplacianT"
        )
     ==
        alphacavg*rhoc*(Uc&g)
      + averagedSh
      + (
            radiationHeat
          ? (radiation->Sh(thermo, he) + fvOptions(rhoc, he))
          : fvOptions(rhoc, he)
        )
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
