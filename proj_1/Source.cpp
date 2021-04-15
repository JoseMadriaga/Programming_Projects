#nclude "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cstdio>
#include <cmath>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

double masses[] = {
	0.00000000000,
	1.00782503223,
	3.0160293201,
	6.0151228874,
	9.012183065,
	10.01293695,
	12.0000000,
	14.00307400443,
	15.99491461957
};

int main()
{
	Molecule mol("/Users/josemarcmadriaga/Documents/prog_proj/proj_1/acetaldehyde.txt", 0);	//calls out the function

	cout << "Number of atoms: " << mol.natom << endl;	// prints number of atoms
	cout << "Input Cartesian coordinates:\n";	//provides a title for the cartesian coordinates
	mol.print_geom();	

	cout << "Interatomic distances (bohr):\n";	//provides a ttle for interatomic distances
	for (int i = 0; i < mol.natom; i++)	// basically goes through unique atom pairs, 1-0, 2-0, 2-1, 3-0 and so forth
		for (int j = 0; j < i; j++)
			printf("%d %d %8.5f\n", i, j, mol.bond(i, j));

	cout << "\nBond angles:\n";		// provides a title for Bond Angle
	for (int i = 0; i < mol.natom; i++) {	// look for unique atomic triplet, 2-1-0, 3-2-1, 3-2-0, 4-3-0, 4-3-1, 4-3-2 but has a condition
		for (int j = 0; j < i; j++) {		// where only those with certain bond distances will have an output to avoid unnecessary bond angles 
			for (int k = 0; k < j; k++) {
				if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0)
					printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i, j, k) * (180.0 / acos(-1.0)));  // (180.0 / acos(-1.0)) radians to degrees
			}
		}
	}

	cout << "\nOut-of-Plane angles:\n";			// provides a title for out of Plane angles
	for (int i = 0; i < mol.natom; i++) {		// look for unique atomic quartet, 4,3,1,0 and has a condition where only tose with certain bond distances will have an output
		for (int k = 0; k < mol.natom; k++) {
			for (int j = 0; j < mol.natom; j++) {
				for (int l = 0; l < j; l++) {
					if (i != j && i != k && i != l && j != k && k != l && mol.bond(i, k) < 4.0 && mol.bond(k, j) < 4.0 && mol.bond(k, l) < 4.0)
						printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.oop(i, j, k, l) * (180.0 / acos(-1.0)));
				}
			}
		}
	}

	cout << "\nTorsional angles:\n\n";		//provides a title for Torsional angles
	for (int i = 0; i < mol.natom; i++) {	// look for atomic quartet that has certain bond distances 
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < j; k++) {
				for (int l = 0; l < k; l++) {
					if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0 && mol.bond(k, l) < 4.0)
						printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.torsion(i, j, k, l) * (180.0 / acos(-1.0)));
				}
			}
		}
	}

	/* find the center of mass (COM) */
	double M = 0.0;		//declaring a total mass to be zero
	for (int i = 0; i < mol.natom; i++) M += masses[(int)mol.zvals[i]];  // adds up each mass corresponding to atomic number

	double xcm = 0.0;	//declaring center of mass x,y,z to zero
	double ycm = 0.0;
	double zcm = 0.0;
	double mi;
	for (int i = 0; i < mol.natom; i++) {
		mi = masses[(int) mol.zvals[i]];
		xcm += mi * mol.geom[i][0];  //sum of mass of each atom * coordinate x
		ycm += mi * mol.geom[i][1];
		zcm += mi * mol.geom[i][2];
	}
	xcm /= M;	// (sum of of each atom * coord x) divided total mass
	ycm /= M;
	zcm /= M;
	printf("\nMolecular center of mass: %12.8f %12.8f %12.8f\n", xcm, ycm, zcm);

	mol.translate(-xcm, -ycm, -zcm); // translating the coordinates to their center of mass

	Matrix I(3, 3);	//constructing a 3 by 3 matrix and declaring them as zeroes
	I.setZero();

	for (int i = 0; i < mol.natom; i++) { //declaring moments of inertia tensor
		mi = masses[(int)mol.zvals[i]];
		I(0, 0) += mi * (mol.geom[i][1] * mol.geom[i][1] + mol.geom[i][2] * mol.geom[i][2]);
		I(1, 1) += mi * (mol.geom[i][0] * mol.geom[i][0] + mol.geom[i][2] * mol.geom[i][2]);
		I(2, 2) += mi * (mol.geom[i][0] * mol.geom[i][0] + mol.geom[i][1] * mol.geom[i][1]);
		I(0, 1) -= mi * mol.geom[i][0] * mol.geom[i][1];
		I(0, 2) -= mi * mol.geom[i][0] * mol.geom[i][2];
		I(1, 2) -= mi * mol.geom[i][1] * mol.geom[i][2];
	}

	I(1, 0) = I(0, 1);	//inertia tensor are symmetrical
	I(2, 0) = I(0, 2);
	I(2, 1) = I(1, 2);

	cout << "\nMoment of inertia tensor (amu bohr^2):\n";	//providing a title for Moment of inertia tensor
	cout << I << endl;

	// Using Eigen to solvefor the eigenvalue and eigenvectors of moments of inertia tensor
	Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
	Matrix evecs = solver.eigenvectors();
	Matrix evals = solver.eigenvalues();

	cout << "\nPrincipal moments of inertia (amu * bohr^2):\n";
	cout << evals << endl;

	double conv = 0.529177249 * 0.529177249;    //bohr to A is 0.529177249
	cout << "\nPrincipal moments of inertia (amu * AA^2):\n";
	cout << evals * conv << endl;

	conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;		// amu to g 1.6605402 x 10^-27 kilograms and A to cm is 10^-8
	cout << "\nPrincipal moments of inertia (g * cm^2):\n";
	cout << evals * conv << endl;

	// classify the rotor 
	if (mol.natom == 2) cout << "\nMolecule is diatomic.\n";
	else if (evals(0) < 1e-4) cout << "\nMolecule is linear.\n";  // linear molecule has a Ia << Ib = Ic where Ia is taken as zero 
	else if ((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))		//spherical top when Ia = Ib = Ic
		cout << "\nMolecule is a spherical top.\n";
	else if ((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) > 1e-4))		// onlate symmetric top Ia= Ib < Ic
		cout << "\nMolecule is an oblate symmetric top.\n";
	else if ((fabs(evals(0) - evals(1)) > 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))		// prolate symmetric top Ia < Ib = Ic 
		cout << "\nMolecule is a prolate symmetric top.\n";
	else cout << "\nMolecule is an asymmetric top.\n";		// different value I's

	// compute the rotational constants 
	// rotational constant expression h/(8*(pi^2)*I)
	double _pi = acos(-1.0);
	conv = 6.6260755E-34 / (8.0 * _pi * _pi);
	conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
	conv *= 1e-6;
	cout << "\nRotational constants (MHz):\n";
	cout << "\tA = " << conv / evals(0) << "\t B = " << conv / evals(1) << "\t C = " << conv / evals(2) << endl;

	// rotational constant expression h/(8*(pi^2)*c*I)
	conv = 6.6260755E-34 / (8.0 * _pi * _pi);
	conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
	conv /= 2.99792458E10; //speed of light in cm/s
	cout << "\nRotational constants (cm-1):\n";
	cout << "\tA = " << conv / evals(0) << "\t B = " << conv / evals(1) << "\t C = " << conv / evals(2) << endl;

	return 0;
}
