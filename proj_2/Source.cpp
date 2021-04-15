#include "molecule.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

int main() {

	Molecule mol("/Users/josemarcmadriaga/Documents/prog_proj/proj_2/benzene.txt", 0);
	mol.print_geom();
	mol.read_hessian("/Users/josemarcmadriaga/Documents/prog_proj/proj_2/benzene_hes.txt");
	mol.mass_weight();
	mol.diag_hessian();
	return 0;
}
