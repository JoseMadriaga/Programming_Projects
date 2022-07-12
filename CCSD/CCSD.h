#include <string> 
#include <cmath> 
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
typedef Eigen::Tensor<double, 4> Tensor4D;

class CCSD
{
public:
	//indices based on basis functions
	int nao;
	int nocc;
	int nvir;
	int nsm;
	int nsm_occ;
	int nsm_vir;
	
	//energy values from SCF, MP2
	double Escf;
	double Emp2;
	double Etot;

	//overlap, kinetic, nuclear, one and two electron integrals; AO = Coefficients
	Matrix overlap;
	Matrix kinetic;
	Matrix nuclear;
	Matrix AO;
	Matrix C_ham;
	std::vector <double> MO_eric;

	//spin orbital two electron integrals, spin orbital Fock Matrix, 
	Tensor4D spi_TMO;
	Matrix spi_Focko;
	
	//t1 and t2 ampltides with corresponding energy denominators 
	Matrix t_ai;
	Tensor4D t_aibj;
	double Dia;
	double Dijab;

	//effective two-particle excitation operators
	Tensor4D tau_eff;
	Tensor4D tau;
	
	//one-particle intermediates
	double Fae;
	double Fmi;
	double Fme;
	
	//two-particle intermediates
	double Wmnij;
	double Wabef;
	double Wmbej;

	//updated t1 and t2 ampltiudes
	Matrix update_tai;
	Tensor4D update_taibj;	

	void read_vals(const char* filename, double &energy);
        void read_txt(const char* filename, Matrix &integrals);
        void read_Matrix(const char* filename, Matrix &core);
        void read_twoe(const char* filename, std::vector <double> &TEI);
        int c_index(int a, int b);
	void conv_spi(std::vector <double> &TEI, Tensor4D &spi_TEI);
	void spi_Fock(Matrix &coeff, Matrix &core, Tensor4D &spi_TEI, Matrix &Fock);
	void guess(Matrix &Fock, Tensor4D &spi_TEI, Matrix &t1_amp, Tensor4D &t2_amp);
	void tau_and_eff(Matrix &t1_amp, Tensor4D &t2_amp, Tensor4D &taue, Tensor4D &ta);		

	double gen_Fae(int a, int e);
	double gen_Fmi(int m, int i);
	double gen_Fme(int m, int e);
	Matrix up_tai();

	double gen_Wmnij(int m, int n, int i, int j);
	double gen_Wabef(int a, int b, int e, int f);
	double gen_Wmbej(int m, int b, int e, int j);
	Tensor4D up_taibj();
	
	double CC_energy(Matrix &t1_amp, Tensor4D &t2_amp);
	//bool converged(Matrix &old_t1amp, Matrix &new_t1amp, Tensor4D &old_t2amp, Tensor4D &new_t2_amp, double &energy_dif);	
	bool converged(double &energy_dif);
};	
