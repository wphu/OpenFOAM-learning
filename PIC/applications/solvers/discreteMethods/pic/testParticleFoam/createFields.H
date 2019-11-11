Info<< nl << "Constructing picCloud " << endl;
IOdictionary picInitialiseDict
(
    IOobject
    (
        "picInitialiseDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary phi_from_file
(
    IOobject
    (
        "electric_potential",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

picCloud pic("pic", mesh, picInitialiseDict, true);

volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);