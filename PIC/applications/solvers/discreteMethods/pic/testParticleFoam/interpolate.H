scalar interpolate(const scalarField& coordinate_z, const scalarField& electric_potential, const scalar& z)
{
    label low = 0;
    label high = coordinate_z.size() - 1;
    label mid = 0;

    if(z <= coordinate_z[low] )
    {
        return electric_potential[low];
    }
    else if(z >= coordinate_z[high])
    {
        /*
        low = high -1;
        double dz_inv = 1.0 / (coordinate_z[low] - coordinate_z[high]);
        return electric_potential[high] * (coordinate_z[low] - z) * dz_inv +
                electric_potential[low] * (z - coordinate_z[high]) * dz_inv;
        */
        return electric_potential[high];
    }

    while(low <= high){
        mid = (low + high) / 2;
        if(coordinate_z[mid] < z){
            low = mid + 1;
        }
        else if(coordinate_z[mid] > z){
            high = mid - 1;
        }
        else {
            return electric_potential[mid];
        }
    }
    // now low is 1 larger than high
    scalar dz_inv = 1.0 / (coordinate_z[low] - coordinate_z[high]);
    return electric_potential[high] * (coordinate_z[low] - z) * dz_inv +
            electric_potential[low] * (z - coordinate_z[high]) * dz_inv;
}
