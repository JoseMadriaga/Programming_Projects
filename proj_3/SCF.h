#include <string>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include <vector>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

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
};
