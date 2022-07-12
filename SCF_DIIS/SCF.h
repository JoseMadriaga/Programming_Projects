#include <string>
#include <iostream> 
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "unsupported/Eigen/CXX11/Tensor"
#include <vector>
#include <numeric>
#include <cmath>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Tensor<double, 3> Tensor3D;
class SCF
{
public:
	Matrix overlap;
	Matrix kinetic; 
	Matrix nuclear; 
	Matrix C_ham;
	Matrix ortho_S;
	Matrix Fock;
	Matrix AO;
	Matrix Int_Dens;
	Matrix old_Dens;
	Matrix new_Dens;
	Matrix temp_Fock;
	double int_E;
	double nuc_e;
	double old_E;
	double new_E;
	double dif_E;
	double Dens_rms;
	int indices;
	int iter;
	int p;
	int q;
	int r;
	int s;
	std::vector <double> twoe;

	void read_txt(const char* filename, Matrix &temp);
	void read_twoe(const char* filename, std::vector <double> &temp);
	void ortho();
	void guess();
	void int_SCF();
	void new_Fock();
	void dif_Dens();
	void iter_SCF();
	void iter_guess();
	void converge();
	void dipole( Matrix &mat, Matrix &mat1, double &temp);
	
	//DIIS
	int num_stored;
	Matrix err;
	Tensor3D err_ten;
	Tensor3D fock_ten;
	Matrix B;
	Matrix ext_Fock;
	void extra_err();
	void stor_errafoc();
	void comp_B(int numb);
};

