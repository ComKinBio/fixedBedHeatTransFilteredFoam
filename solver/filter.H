
IOdictionary CouplingDict
(
    IOobject
    (
        "couplingDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dictionary sourceFilterDict = CouplingDict.subDict("sourceFilter");
const dictionary gasFilterDict = CouplingDict.subDict("gasFilter");

//- make sourceFilter class

Time sourceDiffusionRunTime_
(
    "controlDiffDict",
    args
);

Info
    << "Create source diffusion mesh for time = "
    << sourceDiffusionRunTime_.timeName() << nl << endl;

fvMesh sourceDiffusionMesh_
(
    IOobject
    (
        fvMesh::defaultRegion,
        sourceDiffusionRunTime_.timeName(),
        sourceDiffusionRunTime_,
        IOobject::MUST_READ
    )
);

simpleControl sourceDiffusionSimple_(sourceDiffusionMesh_);

sourceFilter sourceFilterModel
(
    "steelCloud",
    mesh,
    sourceDiffusionRunTime_,
    sourceDiffusionMesh_,
    sourceDiffusionSimple_
);

bool useSourceFilter = readBool(sourceFilterDict.lookup("useSourceFilter"));


//- make gasFilter class

Time gasDiffusionRunTime_
(
    "controlDiffDict",
    args
);

Info
    << "Create gas diffusion mesh for time = "
    << gasDiffusionRunTime_.timeName() << nl << endl;

fvMesh gasDiffusionMesh_
(
    IOobject
    (
        fvMesh::defaultRegion,
        gasDiffusionRunTime_.timeName(),
        gasDiffusionRunTime_,
        IOobject::MUST_READ
    )
);

simpleControl gasDiffusionSimple_(gasDiffusionMesh_);

gasFilter gasFilterModel
(
    mesh,
    gasDiffusionRunTime_,
    gasDiffusionMesh_,
    gasDiffusionSimple_
);

bool useGasFilter = readBool(gasFilterDict.lookup("useGasFilter"));
bool useGasFilterImplicit = readBool(gasFilterDict.lookup("useImplicitLaplacian"));

//- other settings
bool radiationHeat = readBool(CouplingDict.lookup("radiationHeat"));

