Info<< "\nConstructing steelCloud "<< endl;
basicThermoCloud steelParcels
(
    "steelCloud",
    rhoc,
    Uc,
    g,
    slgThermo,
    Ucavg,
    rhoavg,
    muavg,
    alphacavg,
    Tavg,
    kappaavg,
    Cpavg,
    &sourceFilterModel
);

//default value for DEMFlag is 1, if DEM flag is not actived, the solver is normal
//may not need, just for sure, if run from modified cases

steelParcels.DEMFlagSet(1);
