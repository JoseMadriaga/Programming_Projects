#include "SCF.h"
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
	cout << temp << endl;
	txt.close();
	}	
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
		cout << p << '\t' << q << '\t' << r << '\t' << s << '\t' <<  temp[c_index(c_index(p,q),c_index(r,s))] << endl;
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
	for (int i=0; i < indices; i++){
		for (int j=0; j < indices; j++){
			if (abs(ortho_S(i,j)) < 1e-9)
				ortho_S(i,j) = 0;
		}
	}
	cout << ortho_S << endl;
}
void SCF::guess(){
	Fock = ortho_S.transpose() * C_ham * ortho_S;
	        for (int i=0; i < indices; i++){
               		 for (int j=0; j < indices; j++){
                        if (abs(Fock(i,j)) < 1e-9)
                                Fock(i,j) = 0;
                }
        }
	cout << "Fock Matrix" << endl;
	cout << Fock << endl;
	Eigen::SelfAdjointEigenSolver<Matrix> solver(Fock);
	Matrix evecs = solver.eigenvectors();
	Matrix evals = solver.eigenvalues();

	AO = ortho_S*evecs;
	        for (int i=0; i < indices; i++){
                	for (int j=0; j < indices; j++){
                        if (abs(AO(i,j)) < 1e-9)
                                AO(i,j) = 0;
                }
        }
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
	 for (int i=0; i < indices; i++){
                for (int j=0; j < indices; j++){
                        if (abs(Int_Dens(i,j)) < 1e-9)
                                Int_Dens(i,j) = 0;
                }
        }
	cout << "Initial Guess Density" << endl;
	cout << Int_Dens << endl;
			
}
void SCF::int_SCF(){
	int_E = 0;
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
					/*cout << i << j << k << l << twoe[c_index(c_index(i,j),c_index(k,l))] << endl;
					cout << i << j << k << l << twoe[c_index(c_index(i,k),c_index(j,l))] << endl;
					cout << temp_Fock << endl;*/
				}
			}
		}
	}
	/*cout << "New Fock Matrix" << endl;
	cout << temp_Fock << endl;*/
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
        Fock = ortho_S.transpose() * temp_Fock * ortho_S;
                for (int i=0; i < indices; i++){
                         for (int j=0; j < indices; j++){
                        if (abs(Fock(i,j)) < 1e-9)
                                Fock(i,j) = 0;
                }
        }
        /*cout << "Fock Matrix" << endl;
        cout << Fock << endl;*/
        Eigen::SelfAdjointEigenSolver<Matrix> solver(Fock);
        Matrix evecs = solver.eigenvectors();
        Matrix evals = solver.eigenvalues();

        AO = ortho_S*evecs;
                for (int i=0; i < indices; i++){
                        for (int j=0; j < indices; j++){
                        if (abs(AO(i,j)) < 1e-9)
                                AO(i,j) = 0;
                }
        }
        /*cout << "AO matrix" << endl;
        cout << AO << endl;*/

        Int_Dens.resize(indices,indices);
        Int_Dens.setZero();
        for (int i=0; i < indices; i++){
                for (int j=0; j < indices; j++){
                        for (int k=0; k < 5; k++){
                                Int_Dens(i,j) += AO(i,k)*AO(j,k);
                        }
                }
        }
         for (int i=0; i < indices; i++){
                for (int j=0; j < indices; j++){
                        if (abs(Int_Dens(i,j)) < 1e-9)
                                Int_Dens(i,j) = 0;
                }
        }
       /* cout << "Initial Guess Density" << endl;
        cout << Int_Dens << endl;*/
}
void SCF::converge(){
	iter = 1;
	old_E = int_E;
	new_Fock();
	old_Dens = Int_Dens;
	iter_guess();
	new_Dens = Int_Dens;
	iter_SCF();
	new_E = int_E;
	dif_E = abs(new_E - old_E);
	dif_Dens();
	cout << "iterations " << '\t'  << "Total Energy " << '\t' << "Delta E " << '\t' << "RMS " << endl;
        cout << iter << '\t'  << int_E << '\t'  << dif_E << '\t'  << Dens_rms << endl; 
	while ( dif_E > pow(10,-8) && Dens_rms > pow(10,-8) ){
		old_E = int_E;
        	new_Fock();
        	old_Dens = Int_Dens;
        	iter_guess();
        	new_Dens = Int_Dens;
        	iter_SCF();
        	new_E = int_E;
        	dif_E = abs(new_E - old_E);
		dif_Dens();
		iter += 1;
		cout << "iterations" << '\t' << "Total Energy" << '\t' << "Delta E" << '\t' << "RMS" << endl;
        	cout << iter << '\t' << int_E << '\t' << dif_E << '\t' << Dens_rms << endl; 
	}
	cout << "Final Energy" << endl;
	printf("%8.12f/n", int_E);
}
