#include "SCF.h"
#include "molecule.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
using namespace std;

void SCF::read_txt(const char* filename,Matrix &temp) {
	ifstream txt(filename);
	assert(txt.good());
	double num1 = 0;
	double num2 = 0; 
	temp.resize(indices, indices);
	temp.setZero();
	while(txt)
	{
		for (int i = 0; i < indices; i++) {
			for (int j = 0; j < i+1; j++) {
				txt >> num1 >> num2 >> temp(i,j);
				temp(j,i)=temp(i,j);
			}
		}
	}	
	cout << temp << endl;
	txt.close();
}

int c_index(int a, int b){
	int val=0;
	if (a>b)
	{
		val = a*(a+1)/2 + b;
	}
	else
	{
		val = b*(b+1)/2 + a;
	}
	return val; 
}
void SCF::read_twoe(const char* filename, vector <double>  &temp) {
	ifstream txt(filename);
	assert(txt.good());
	temp.resize(26106);
	while(txt)
	{
		txt >> p >> q >> r >> s >> temp[c_index(c_index(p,q),c_index(r,s))];
		cout << c_index(c_index(p,q),c_index(r,s))  << '\t' << p << '\t' << q << '\t' << r << '\t' << s << '\t' <<  temp[c_index(c_index(p,q),c_index(r,s))] << endl;
	}
	txt.close();
}

void SCF::ortho() {
	Eigen::SelfAdjointEigenSolver<Matrix> solver(overlap);
	Matrix evecs = solver.eigenvectors();
	Matrix evals = solver.eigenvalues();

	Matrix evals_mat;

	evals_mat.resize(indices,indices);
	evals_mat.setZero();
	for (int i = 0; i < indices; i++){
		evals_mat(i,i) = 1/sqrt(evals(i));
	}
	ortho_S = evecs*evals_mat*evecs.transpose();
	cout << ortho_S << endl;
}
void SCF::guess(){
	Fock = ortho_S.transpose() * C_ham * ortho_S;
	cout << "Fock Matrix" << endl;
	cout << Fock << endl;
	Eigen::SelfAdjointEigenSolver<Matrix> solver(Fock);
	Matrix evecs = solver.eigenvectors();
	Matrix evals = solver.eigenvalues();

	AO = ortho_S*evecs;
	cout << "AO matrix" << endl;
	cout << AO << endl;
	
	Int_Dens.resize(indices,indices);
	Int_Dens.setZero();
	for (int i=0; i < indices; i++){
		for (int j=0; j < indices; j++){
			for (int k=0; k < 5; k++){
				Int_Dens(i,j) += AO(i,k)*AO(j,k);
			}
		}
	}
	cout << "Initial Guess Density" << endl;
	cout << Int_Dens << endl;
			
}
void SCF::int_SCF(){
	int_E = 0;
	//temp_Fock.setZero(indices,indices);
	//temp_Fock = C_ham;
	for (int i = 0; i < indices; i++){
		for (int j = 0; j < indices; j++){
			int_E += Int_Dens(i,j)*(C_ham(i,j)*2);
		}
	}
	int_E += nuc_e;
	cout << int_E << endl;
}
void SCF::new_Fock(){
	temp_Fock.resize(indices,indices);
        temp_Fock.setZero();
	for (int i = 1; i < indices+1; i++){
		for (int j=1; j < indices+1; j++){
			temp_Fock(i-1,j-1) += C_ham(i-1,j-1);
			for (int k=1; k < indices+1; k++){
				for (int l=1; l < indices+1; l++){
					temp_Fock(i-1,j-1) += Int_Dens(k-1,l-1)*(2*twoe[c_index(c_index(i,j),c_index(k,l))] - twoe[c_index(c_index(i,k),c_index(j,l))]);
				}
			}
		}
	}
}
void SCF::dif_Dens(){
	Dens_rms = 0;
	for (int i =0; i < indices; i++){
		for (int j=0; j < indices; j++){
			Dens_rms+= pow(new_Dens(i,j)-old_Dens(i,j),2.0);
		}
	}
	Dens_rms = sqrt(Dens_rms);
}

void SCF::iter_SCF(){
        int_E = 0;
        for (int i = 0; i < indices; i++){
                for (int j = 0; j < indices; j++){
                        int_E += Int_Dens(i,j)*(C_ham(i,j)+temp_Fock(i,j));
                }
        }
        int_E += nuc_e;
        /*cout << int_E << endl;*/
}
void SCF::iter_guess(){
	if (iter > 1)
		temp_Fock = ext_Fock;
        Fock = ortho_S.transpose() * temp_Fock * ortho_S;    
        Eigen::SelfAdjointEigenSolver<Matrix> solver(Fock);
        Matrix evecs = solver.eigenvectors();
        Matrix evals = solver.eigenvalues();
        AO = ortho_S*evecs;
        Int_Dens.resize(indices,indices);
        Int_Dens.setZero();
        for (int i=0; i < indices; i++){
                for (int j=0; j < indices; j++){
                        for (int k=0; k < 5; k++){
                                Int_Dens(i,j) += AO(i,k)*AO(j,k);
                        }
                }
        }
}

// DIIS procedure 
/* First, compute the error matrices, e_i, and store it, 
 Second,build the B matrix which is solve using Householder
Third, generate the new fock matrix F' = ciFi
*/
//calculates the error with each iterative 
void SCF::extra_err(){
	err.setZero(indices,indices);
	err = temp_Fock*Int_Dens*overlap - overlap*Int_Dens*temp_Fock;
	//err = ortho_S.transpose() * err * ortho_S;
	//cout << ortho_S.transpose() << endl;
}

//stores the error matrices and the fock matrices	
void SCF::stor_errafoc(){
	for (int j=0; j < indices; j++){
		for (int k = 0; k < indices; k++){
				int i = (iter-1) % 6;			
				extra_err();
				err_ten(i,j,k) = err(j,k);
				fock_ten(i,j,k) = temp_Fock(j,k);		
				}
		}
	cout << "error matrix" << endl;
	cout << (iter-1) %6 << '\t' << err << endl;
	cout << "fock matrix" << endl;
	cout << (iter-1)%6 << '\t' << temp_Fock << endl;
	
/*	array<long,3> offset = {0,0,0};         //Starting point
        array<long,3> extent = {1 ,7,7};       //Finish point
        array<long,2> shape2 = {7,7};         //Shape of desired rank-2 tensor (matrix)
        cout << fock_ten.slice(offset,extent).reshape(shape2) << endl;
*/

}

//build the B matrix and solve for C

void SCF::comp_B(int numb){
		//cout << "here" << endl;
		if (numb > 5)
			numb = 6;
		B.resize(numb+1,numb+1);
		for (int i=0; i < numb; i++)
			for (int j=0; j < numb; j++){
					B(i,j) = 0;
					for (int k=0; k < indices; k++)
						for (int l=0; l < indices; l++){
							B(i,j) += err_ten(i,k,l)*err_ten(j,k,l);
							//cout << i << '\t' << j << '\t' << k << '\t' << l << endl;
                         				//cout << B(i,j) << '\t' << err_ten(i,k,l) << '\t' << err_ten(j,k,l) << endl;
						}
					B(i,numb) = -1;
					B(numb, j) = -1;
					}
		B(numb,numb) = 0; 
		cout << "B" << endl;
		cout << B << endl;
		Eigen::VectorXd C;
		C.setZero(numb+1);
		C(numb) = -1;
		cout << "C" << endl;
		cout << C << endl;
                Eigen::HouseholderQR<Matrix> mat(B);
		Eigen::VectorXd ans;
                ans = mat.solve(C);
            	cout << ans << endl;

		ext_Fock.setZero(indices,indices);
		for (int i=0; i < numb; i++)
			for (int j=0; j < indices; j++)
				for (int k=0; k < indices; k++)
					ext_Fock(j,k) += ans(i)*fock_ten(i,j,k);
}	

				
void SCF::converge(){
	iter = 1;
	err_ten.resize(6,indices, indices);
        fock_ten.resize(6,indices, indices);	
	old_E = int_E;
	new_Fock();
	cout << "initial Fock Matrix" << endl;
	cout << temp_Fock << endl;
	old_Dens = Int_Dens;
	stor_errafoc();
	iter_guess();
	new_Dens = Int_Dens;
	iter_SCF();
	new_E = int_E;
	dif_E = abs(new_E - old_E);
	dif_Dens();
	cout << "iterations " << '\t'  << "Total Energy " << '\t' << "Delta E " << '\t' << "RMS " << endl;
        cout << iter << '\t'  << int_E << '\t'  << dif_E << '\t'  << Dens_rms << endl; 
	printf("%8.12f \n", int_E);
	while ( dif_E > pow(10,-13) && Dens_rms > pow(10,-13) ){
		iter += 1;
		old_E = int_E;
        	new_Fock();
        	old_Dens = Int_Dens;
        	stor_errafoc();
		if (iter > 1)
			comp_B(iter);	
		iter_guess();
        	new_Dens = Int_Dens;
        	iter_SCF();
        	new_E = int_E;
        	dif_E = abs(new_E - old_E);
		dif_Dens();
		cout << "iterations" << '\t' << "Total Energy" << '\t' << "Delta E" << '\t' << "RMS" << endl;
        	cout << iter << '\t' << int_E << '\t' << dif_E << '\t' << Dens_rms << endl; 
		printf(" %8.12f \n", int_E); 
	}
	cout << "Final Energy" << endl;
	printf("%8.12f", int_E);

}
 
