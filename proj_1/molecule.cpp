#include "molecule.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

// prints out the coordinates with corresponding atomic number, Z
void Molecule::print_geom()
{
    for (int i = 0; i < natom; i++)
        printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

//translate the coordinates
void Molecule::translate(double x, double y, double z)
{
    for (int i = 0; i < natom; i++) {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

//looks for the file, gets the charge, declare variables 
Molecule::Molecule(const char* filename, int q) {
    charge = q;

    // open filename
    std::ifstream is(filename);
    assert(is.good());

    // read the number of atoms from filename
    is >> natom;

    // allocate space
    zvals = new int[natom];
    geom = new double* [natom];
    for (int i = 0; i < natom; i++)
        geom[i] = new double[3];

    for (unsigned int i = 0; i < natom; i++)
        is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

    is.close();
}

// bond distance formula
double Molecule::bond(int a, int b)
{
    return sqrt((geom[a][0] - geom[b][0]) * (geom[a][0] - geom[b][0])
        + (geom[a][1] - geom[b][1]) * (geom[a][1] - geom[b][1])
        + (geom[a][2] - geom[b][2]) * (geom[a][2] - geom[b][2]));
}
//converting bond distance to unit vectors
double Molecule::unit(int cart, int a, int b)
{
    return -(geom[a][cart] - geom[b][cart]) / bond(a, b);
}

// Returns the angle between atoms a, b, and c in radians
double Molecule::angle(int a, int b, int c)
{
    return acos(unit(0, b, a) * unit(0, b, c) + unit(1, b, a) * unit(1, b, c) + unit(2, b, a) * unit(2, b, c));
}

//returns out of plane 
// expression is asin((ekj x ekl) * eki)/sin(theta jkl))
double Molecule::oop(int a, int b, int c, int d)
{
    // cross product of unit vector cb and cd
    // normal vector of plane bcd
    double ebcd_x = (unit(1, c, b) * unit(2, c, d) - unit(2, c, b) * unit(1, c, d));
    double ebcd_y = (unit(2, c, b) * unit(0, c, d) - unit(0, c, b) * unit(2, c, d));
    double ebcd_z = (unit(0, c, b) * unit(1, c, d) - unit(1, c, b) * unit(0, c, d));

    //dot product of normal vector of place bcd * unit vector ca
    double exx = ebcd_x * unit(0, c, a);
    double eyy = ebcd_y * unit(1, c, a);
    double ezz = ebcd_z * unit(2, c, a);

    //sum of dot product x,y,z divided by the sin of angle(j,k,l)
    double theta = (exx + eyy + ezz) / sin(angle(b, c, d));

    //numerical precision if value is out of bond 
    if (theta < -1.0) theta = asin(-1.0);
    else if (theta > 1.0) theta = asin(1.0);
    else theta = asin(theta);

    return theta;
}

// Computes the angle between planes a-b-c and b-c-d
// expression acos(((eij x ejk) * (ejk x ekl))/(sin(theta ijk) * sin(theta jkl)))
double Molecule::torsion(int a, int b, int c, int d)
{
    //normal vector of plane abc
    double eabc_x = (unit(1, b, a) * unit(2, b, c) - unit(2, b, a) * unit(1, b, c));
    double eabc_y = (unit(2, b, a) * unit(0, b, c) - unit(0, b, a) * unit(2, b, c));
    double eabc_z = (unit(0, b, a) * unit(1, b, c) - unit(1, b, a) * unit(0, b, c));
    
    //normal vector of plane bcd
    double ebcd_x = (unit(1, c, b) * unit(2, c, d) - unit(2, c, b) * unit(1, c, d));
    double ebcd_y = (unit(2, c, b) * unit(0, c, d) - unit(0, c, b) * unit(2, c, d));
    double ebcd_z = (unit(0, c, b) * unit(1, c, d) - unit(1, c, b) * unit(0, c, d));

    //dot product of the normal vectors
    double exx = eabc_x * ebcd_x;
    double eyy = eabc_y * ebcd_y;
    double ezz = eabc_z * ebcd_z;

    //sum of x,y,z and normalized
    double tau = (exx + eyy + ezz) / (sin(angle(a, b, c)) * sin(angle(b, c, d)));
    
    //numerical precision for out of bound
    if (tau < -1.0) tau = acos(-1.0);
    else if (tau > 1.0) tau = acos(1.0);
    else tau = acos(tau);

    // Compute the sign of the torsion 
    // expression is cross product normal vector of normal vectors abc and bcd
    double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
    double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
    double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
    double norm = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
    cross_x /= norm;
    cross_y /= norm;
    cross_z /= norm;
    double sign = 1.0;
    // cross product times unit vector bc
    double dot = cross_x * unit(0, b, c) + cross_y * unit(1, b, c) + cross_z * unit(2, b, c);
    if (dot < 0.0) sign = -1.0;

    return tau * sign;
}

// deletes the variables
Molecule::~Molecule() {
    delete[] zvals;
    for (int i = 0; i < natom; i++)
        delete[] geom[i];
    delete[] geom;
}
