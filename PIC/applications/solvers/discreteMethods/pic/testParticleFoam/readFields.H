scalarField coordinate_z(phi_from_file.lookup("coordinate_z"));
scalarField electric_potential(phi_from_file.lookup("electric_potential"));
Info<<coordinate_z<<"  "<<electric_potential<<endl;


scalarField& phi_internal = phi.primitiveFieldRef();
const volVectorField& C = mesh.C();
forAll(phi, celli)
{
    phi_internal[celli] = interpolate(coordinate_z, electric_potential, C[celli].component(2));
}
phi.correctBoundaryConditions();

E=-fvc::grad(phi);


