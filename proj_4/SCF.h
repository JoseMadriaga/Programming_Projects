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
	Matrix MO_coef;
	Matrix Int_Dens;
	Matrix old_Dens;
	Matrix new_Dens;
	Matrix temp_Fock;
	Matrix MO_Fock;
	Matrix eMO_Fock;
	Matrix mux;
	Matrix muy;
	Matrix muz;
	double dip_mom;
	double int_E;
	double nuc_e;
	double old_E;
	double new_E;
	double dif_E;
	double Dens_rms;
	double E_MP2;
	double smarE_MP2;
	int indices;
	int iter;
	int p;
	int q;
	int r;
	int s;
	int ijkl;
	std::vector <double> twoe;
	std::vector <double> twoe_MO;
	std::vector <double> smartwoe_MO;	
	
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
	void MO();
	void eigen_MO();
	void conv_twoe();
	void smarconv_twoe();
	void MP2_energy();
	void smarMP2_energy();
};
