#include "SCF.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using namespace std;

SCF fun;

int main() {

	ifstream nuc("/Users/josemarcmadriaga/Documents/prog_proj/proj_3/enuc.txt");
	assert(nuc.good());

	nuc >> fun.nuc_e;
	cout << fun.nuc_e << '\n';

	nuc.close();

	fun.indices = 7;
	cout << "overlap matrix" << '\n';
	fun.read_txt("/Users/josemarcmadriaga/Documents/prog_proj/proj_3/s.txt", fun.overlap);
	cout << "kinetic matrix" << '\n';
	fun.read_txt("/Users/josemarcmadriaga/Documents/prog_proj/proj_3/t.txt",fun.kinetic);
	cout << "nuclear-attraction matrix" << '\n';
	fun.read_txt("/Users/josemarcmadriaga/Documents/prog_proj/proj_3/v.txt", fun.nuclear);

	fun.C_ham.resize(fun.indices, fun.indices);
	fun.C_ham.setZero();
	fun.C_ham = fun.kinetic + fun.nuclear;
	cout << "Core Hamiltonian" << '\n';
	cout << fun.C_ham << '\n';

	cout << "Two electron integral" << '\n';
	fun.read_twoe("/Users/josemarcmadriaga/Documents/prog_proj/proj_3/eric.dat", fun.twoe);

	cout << "Orthogonalization matrix" << '\n';
	fun.ortho();

	fun.guess();

	fun.int_SCF();

	fun.converge();	
		
	fun.MO();
	cout << " Two electron integral in MO" << endl; 
	//fun.conv_twoe();
	fun.smarconv_twoe();
	cout << "SCF_energy" << endl;
	printf(" %2.12f\n ", fun.int_E);
	cout << "MP2 energy" << endl;
	//fun.MP2_energy();
	fun.smarMP2_energy();
	cout << "Total Energy" << endl;
	//printf(" %2.12f\n ", fun.int_E + fun.E_MP2);
	printf(" %2.12f\n ", fun.int_E + fun.smarE_MP2);
	return 0;
};
