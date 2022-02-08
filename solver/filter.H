
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
const dictionary TFMDict = CouplingDict.subDict("TFM");
const dictionary DEMFlagDict = CouplingDict.subDict("DEMFlag");

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


//- make TFM class

bool useTFMGas = readBool(TFMDict.lookup("useTFMGas"));
bool useTFMSource = readBool(TFMDict.lookup("useTFMSource"));

vector gridOrigin;
TFMDict.lookup("gridOrigin") >> gridOrigin;
vector gridVertex;
TFMDict.lookup("gridVertex") >> gridVertex;
scalar gridSize = readScalar(TFMDict.lookup("gridSize"));

coarserGrid TFM
(
    mesh,
    gridOrigin,
    gridVertex,
    gridSize
);


// //- DEMFlag settings
// bool useDEMFlag = readBool(DEMFlagDict.lookup("useDEMFlag"));
// //in every fluid time step, DEM invoked from 0 to Max step (DEMStepMax)
// label DEMStep = 0;
// label DEMStepMax = readLabel(DEMFlagDict.lookup("DEMStepMax")); //requires input
// 
// //DEMInter is the fluid time step to reset DEMFlag again
// scalar DEMInter = readScalar(DEMFlagDict.lookup("DEMInter"));//requires input
// scalar DEMNextTime = 0.0;


//- other settings
bool radiationHeat = readBool(CouplingDict.lookup("radiationHeat"));

