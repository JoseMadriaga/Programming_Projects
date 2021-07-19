#include "CCSD.h"

using namespace std;

//reads a file and stores an energy value

void CCSD::read_vals(const char* filename, double &energy){

	ifstream vals(filename);
	assert(vals.good());

	vals >> energy;

	vals.close();
}

//reads a file containing three columns in order of AO basis functions, AO basis functions, and integral values and stores them in a matrix

void CCSD::read_txt(const char* filename, Matrix &integrals) {
	ifstream txt(filename);
	assert(txt.good());
	double num1=0;
	double num2=0;
	integrals.setZero(nao,nao);
	while(txt)
	{
		for(int i=0; i <nao; i++)
			for (int j=0; j<i+1; j++){
				txt >> num1 >> num2 >> integrals(i,j);
				integrals(j,i)=integrals(i,j);
			}
	
	txt.close();
	}
}

//reads a file containing the core Hamiltonian matrix and stores them in a matrix

void CCSD::read_Matrix(const char* filename, Matrix &core){
	ifstream cor(filename);
	assert(cor.good());
	double num1 = 0;
	double num2 = 0;
	core.setZero(nao,nao);
	while(cor)
	{
		for (int i=0; i< nao; i++)
			for (int j=0; j < nao; j++){
				cor >> num1 >> num2 >> core(i,j);
			}
	cor.close();
	}
}

//transforming a four dimensional array into a single array using compound indices 

int CCSD::c_index(int a,int b){
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

//reads a file containing an array of two electron integral 

void CCSD::read_twoe(const char* filename, vector <double> &TEI){
	ifstream MO(filename);
	assert(MO.good());
	TEI.resize(pow(nao,4));
	int p,q,r,s;
	while(MO)
	{
		MO >> p >> q >> r >> s >> TEI[c_index(c_index(p,q),c_index(r,s))];
	}
	MO.close();
}
	
//Converts spatial MO two electron integral to spin MO antisymmetrized two electron integral

void CCSD::conv_spi(vector <double> &TEI, Tensor4D &spi_TEI){
        spi_TEI	= Tensor4D (nsm,nsm,nsm,nsm);
	int pr, qs, prqs, ps, qr, psqr;
        double val1, val2;
        for (int p=0; p < nsm; p++)
                for(int q=0; q < nsm; q++)
                        for (int r=0; r < nsm; r++)
                                for (int s=0; s < nsm; s++){
                                        pr = c_index(p/2,r/2);
                                        qs = c_index(q/2,s/2);
                                        prqs = c_index(pr,qs);
                                        val1 = TEI[prqs] * (p%2 == r%2) * (q%2 == s%2);
                                        ps = c_index(p/2,s/2);
                                        qr = c_index(q/2,r/2);
                                        psqr = c_index(ps,qr);
                                        val2 = TEI[psqr] * (p%2 == s%2) * (q%2 == r%2);
                                        spi_TEI(p,q,r,s) = val1 - val2;
					//if (abs(spi_TEI(p,q,r,s)) < 1e-12){
						//spi_TEI(p,q,r,s) = 0;
					//}
                                        //cout << p << '\t' << q << '\t' << r << '\t' << s << '\t' << spi_TEI(p,q,r,s) << '\t' << val1 << '\t' << val2 << endl;
                                }
}

//Converts AO Core Hamiltonian to spin MO and is added to the spin MO Fock Matrix 
//Adds the spin MO antisymmetrized two electron integral to spin MO fock Matrix

void CCSD::spi_Fock(Matrix &coeff, Matrix &core, Tensor4D &spi_TEI, Matrix &Fock){
	Fock.setZero(nsm,nsm);
        core = coeff.transpose() * core * coeff;
        for ( int p = 0; p < nsm; p++)
                for ( int q = 0; q < nsm; q++){
                	Fock(p,q) = core(p/2,q/2) * (p%2 == q%2);
       	        }
        //cout << "spin orbital Fock Matrix" << endl;
        for (int p=0; p < nsm; p++)
                for ( int q=0; q < nsm; q++)
                        for (int m=0; m < nsm_occ; m++){
                                Fock(p,q) += spi_TEI(p,m,q,m);
				if (abs(Fock(p,q)) < 1e-5){
					Fock(p,q) = 0;
                        
				}
			}

        //cout << Fock << endl;
}

//diagonalizes the spin MO Fock Matrix, creates the initial guess for t1 and t2 amplitudes, and check if the Emp2 generated from this initial guess is the same as the Emp2 from MP2 project

void CCSD::guess(Matrix &Fock, Tensor4D &spi_TEI, Matrix &t1_amp, Tensor4D &t2_amp){
        Eigen::SelfAdjointEigenSolver<Matrix> solver(Fock);
        Matrix evals = solver.eigenvalues();
        double Emp2_int;
	t1_amp.setZero(nsm_occ, nsm_vir);
        t2_amp = Tensor4D (nsm_occ, nsm_occ, nsm_vir, nsm_vir);
        for ( int i = 0; i < nsm_occ; i++)
                for ( int j=0; j < nsm_occ; j++)
                        for (int a=0; a < nsm_vir; a++)
                                for (int b=0; b < nsm_vir; b++){
                                        Dijab = evals(i) + evals(j) - evals(a+nsm_occ) - evals(b+nsm_occ);
                                        t2_amp(i,j,a,b) = spi_TEI(i,j,a+nsm_occ,b+nsm_occ)/Dijab;
                                        //Emp2_int += spi_TEI(i,j,a+nsm_occ,b+nsm_occ)*t2_amp(i,j,a,b);
                                        }
	//cout << Emp2_int/4 << endl;
}

//creates the effective tau and tau 

void CCSD::tau_and_eff(Matrix &t1_amp, Tensor4D &t2_amp, Tensor4D &taue, Tensor4D &ta){
	//cout << "taue" << endl;
	taue = Tensor4D (nsm_occ, nsm_occ, nsm_vir, nsm_vir);
	taue = t2_amp;
	Tensor4D taue_2(nsm_occ, nsm_occ, nsm_vir, nsm_vir);
        for (int i=0; i < nsm_occ; i++)
                for (int j=0; j < nsm_occ; j++)
                        for (int a=0; a < nsm_vir; a++)
                                for (int b=0; b < nsm_vir; b++){
                                        taue_2(i,j,a,b) += t1_amp(i,a)*t1_amp(j,b) - t1_amp(i,b)*t1_amp(j,a);
                                        taue(i,j,a,b) += 0.5 * taue_2(i,j,a,b);
                                        //cout << taue(i,j,a,b) << endl;
                                }
	//cout << "ta" << endl;
	ta = Tensor4D (nsm_occ, nsm_occ, nsm_vir, nsm_vir);
	tau = t2_amp;
        for (int i=0; i < nsm_occ; i++)
                for (int j=0; j < nsm_occ; j++)
                        for (int a=0; a< nsm_vir; a++)
                                for (int b=0; b<nsm_vir ; b++){
                                        ta(i,j,a,b) += t1_amp(i,a)*t1_amp(j,b) - t1_amp(i,b)*t1_amp(j,a);
					//cout << ta(i,j,a,b) << endl;
                                }
} 

//creates one-particle intermediates Fae, Fmi, Fme
//needs the parameters mentioned in lowercase of after the F of each term
double CCSD::gen_Fae(int a, int e){	
	double Fae_2 = 0;
	double Fae_3 = 0;
	double Fae_4 = 0;
        for (int m=0; m < nsm_occ; m++){
        	Fae_2 += spi_Focko(m,e+nsm_occ)*t_ai(m,a);
               	for (int f=0; f < nsm_vir; f++){
                	Fae_3 += t_ai(m,f)*spi_TMO(m,a+nsm_occ, f+nsm_occ, e+nsm_occ);
                        	for (int n=0; n < nsm_occ; n++){
                                	Fae_4 += tau_eff(m,n,a,f)*spi_TMO(m,n,e+nsm_occ, f+nsm_occ);
                                }
                        }
                }
        Fae = (1 - a==e )*spi_Focko(a+nsm_occ,e+nsm_occ) - Fae_2/2 + Fae_3 - Fae_4/2;
       	return Fae;
}

double CCSD::gen_Fmi(int m, int i){
        double Fmi = 0;
        double Fmi_2 = 0;
        double Fmi_3 = 0;
        double Fmi_4 = 0;
        for (int e=0; e < nsm_vir; e++){
       		Fmi_2 += t_ai(i,e)*spi_Focko(m,e+nsm_occ);
                for (int n=0; n < nsm_occ; n++){
                	Fmi_3 += t_ai(n,e)*spi_TMO(m,n,i,e+nsm_occ);
                        for (int f=0; f<nsm_vir ; f++){
                        	Fmi_4 += tau_eff(i,n,e,f)*spi_TMO(m,n,e+nsm_occ, f+nsm_occ);
                        }
                }
	}
       	Fmi = (1 - m==i )*spi_Focko(m,i) + Fmi_2/2 + Fmi_3 + Fmi_4/2;
	return Fmi;
}
double CCSD::gen_Fme(int m, int e){
       	double Fme = 0;
	double Fme_2 = 0;
	for (int n=0; n < nsm_occ; n++)
        	for (int f=0; f< nsm_vir; f++){
                        Fme_2 += t_ai(n,f)*spi_TMO(m,n,e+nsm_occ,f+nsm_occ);
                }
        Fme = spi_Focko(m,e+nsm_occ) + Fme_2;
        return Fme;
}

//Update t1 amplitudes
Matrix CCSD::up_tai(){
	Eigen::SelfAdjointEigenSolver<Matrix> solver(spi_Focko);
        Matrix evals = solver.eigenvalues();	
	double t1_1, t1_2, t1_3, t1_4, t1_5, t1_6, t1_7;
	update_tai.setZero(nsm_occ,nsm_vir);
        for (int i=0; i< nsm_occ; i++)
                for (int a=0; a<nsm_vir; a++){
			t1_2 = 0;
			for (int e=0; e< nsm_vir; e++){
		
                                t1_2 += t_ai(i,e)*gen_Fae(a,e);
                	}
			t1_3 = 0;
			t1_4 = 0;
			t1_6 = 0;
			t1_7 = 0;
		        for (int m=0; m< nsm_occ; m++){
        
	                	t1_3 += t_ai(m,a)*gen_Fmi(m,i);
        
	                	for (int e=0; e < nsm_vir; e++){
 			        	
					t1_4 += t_aibj(i,m,a,e)*gen_Fme(m,e);
                                        
					for (int f=0; f< nsm_vir; f++){
                                               	
						t1_6 += t_aibj(i,m,e,f)*spi_TMO(m,a+nsm_occ,e+nsm_occ,f+nsm_occ);
						
                       				if (i == 4 && a == 2){
							//cout << "t1_6(4,2)" << '\t' <<  i << '\t' << a << '\t' << t_aibj(i,m,e,f) << '\t' << spi_TMO(m,a+nsm_occ,e+nsm_occ, f+nsm_occ) << endl;
		                 		}
					}        
						
					for (int n=0; n<nsm_occ; n++){
                                                        
						t1_7 += t_aibj(m,n,a,e)*spi_TMO(n,m,e+nsm_occ,i);
						if (i ==4 && a ==2){
							//cout << "t1_7(4,2)" << '\t' << i << '\t' << a << '\t' << t_aibj(m,n,a,e) << '\t' << spi_TMO(n,m,e+nsm_occ,i) << endl;
                                        	}
					}

				}
			}
			t1_5 = 0;
			for (int n=0; n< nsm_occ; n++)
				for (int f=0; f< nsm_vir; f++){
					
					t1_5 += t_ai(n,f)*spi_TMO(n,a+nsm_vir,i,f+nsm_occ);
			
                                        }
                        Dia = evals(i) - evals(a+nsm_occ);
			t1_1 += spi_Focko(i,a+nsm_occ); 
			update_tai(i,a) = t1_1 + t1_2 - t1_3 + t1_4 - t1_5 - t1_6/2 - t1_7/2;
                       	//cout << i << '\t' << a << '\t' << -t1_6/2 << '\t' << - t1_6/2 - t1_7/2 << '\t' << (-t1_6/2 - t1_7/2)/Dia << endl; 
			update_tai(i,a) /= Dia;
			if (abs(update_tai(i,a)) < 1e-10){
				update_tai(i,a) = 0;
			}
                }
        return update_tai;
}

//creates two-particle intermediates Wmnij, Wabef, Wmbej

double CCSD::gen_Wmnij(int m, int n, int i, int j){
	double Wmnij_2 = 0;
	double Wmnij_3 = 0;
        for (int e=0; e<nsm_vir; e++){
        	Wmnij_2 += t_ai(j,e) * spi_TMO(m,n,i,e+nsm_occ) - t_ai(i,e) * spi_TMO(m,n,j,e+nsm_occ);
		for (int f=0; f<nsm_vir; f++){
	        	Wmnij_3 += tau(i,j,e,f)*spi_TMO(m,n,e+nsm_occ,f+nsm_occ);
                }
	}
        Wmnij = spi_TMO(m,n,i,j) + Wmnij_2 + Wmnij_3/4;   
	return Wmnij;
}

double CCSD::gen_Wabef(int a, int b, int e, int f){
        double Wabef_2 = 0;
        double Wabef_3 = 0;
	for (int m=0; m < nsm_occ; m++){
        	Wabef_2 += t_ai(m,b)*spi_TMO(a+nsm_occ,m,e+nsm_occ,f+nsm_occ);
		Wabef_2 -= t_ai(m,a)*spi_TMO(b+nsm_occ,m,e+nsm_occ,f+nsm_occ);
		for (int n=0; n< nsm_occ; n++){
                	Wabef_3 += tau(m,n,a,b)*spi_TMO(m,n,e+nsm_occ,f+nsm_occ);
                }
	}
	Wabef = spi_TMO(a+nsm_occ,b+nsm_occ,e+nsm_occ,f+nsm_occ) - Wabef_2 + Wabef_3/4;
	return Wabef;
}


double CCSD::gen_Wmbej(int m, int b, int e, int j){
        double Wmbej_2 = 0;
        double Wmbej_3 = 0;
        double Wmbej_4 = 0;
	for (int f=0; f<nsm_vir; f++){
		Wmbej_2 += t_ai(j,f)*spi_TMO(m,b+nsm_occ,e+nsm_occ,f+nsm_occ);
        }
	for (int n=0; n < nsm_occ; n++){
                Wmbej_3 += t_ai(n,b)*spi_TMO(m,n,e+nsm_occ,j);
		for (int f=0; f < nsm_vir; f++){
                        Wmbej_4 += (0.5 * t_aibj(j,n,f,b) + t_ai(j,f)*t_ai(n,b))*spi_TMO(m,n,e+nsm_occ,f+nsm_occ);
                }
       	}
       	Wmbej = spi_TMO(m,b+nsm_occ,e+nsm_occ,j) + Wmbej_2 - Wmbej_3 - Wmbej_4;
	return Wmbej;
}

// Updates t2 ampltiudes
Tensor4D CCSD::up_taibj(){
	Eigen::SelfAdjointEigenSolver<Matrix> solver(spi_Focko);
        Matrix evals = solver.eigenvalues();
	double t2_1, t2_2, t2_3, t2_4, t2_5, t2_6, t2_7, t2_8;
        update_taibj = Tensor4D (nsm_occ,nsm_occ,nsm_vir,nsm_vir);
        for (int i=0; i< nsm_occ; i++)
                for (int j=0; j< nsm_occ; j++)
                        for (int a=0; a<nsm_vir; a++)
                                for (int b=0; b<nsm_vir; b++){
					t2_2 = 0;
					t2_5 = 0;
					t2_7 = 0;
					for (int e=0; e< nsm_vir;e++){
						double t2_2a = 0;
						double t2_2b = 0;
						for (int m=0; m < nsm_occ; m++){
							t2_2a += t_ai(m,b)*gen_Fme(m,e);
							t2_2b += t_ai(m,a)*gen_Fme(m,e);
	 					}
						t2_2 += t_aibj(i,j,a,e)*(gen_Fae(b,e) - 0.5*t2_2a);
						t2_2 -= t_aibj(i,j,b,e)*(gen_Fae(a,e) - 0.5*t2_2b);
						//cout << "em" << '\t' << e << endl;
						for (int f=0; f< nsm_vir; f++){
							t2_5 += tau(i,j,e,f)*gen_Wabef(a,b,e,f);
							//cout << "ef" << '\t' << e << '\t' << f << endl;
						}	
						t2_7 += t_ai(i,e)*spi_TMO(a+nsm_occ,b+nsm_occ,e+nsm_occ,j);
						t2_7 -= t_ai(j,e)*spi_TMO(a+nsm_occ,b+nsm_occ,e+nsm_occ,i);
						//cout << "e" << '\t' << e  << endl;
					}
					
					t2_3 = 0;
					t2_4 = 0;
					t2_6 = 0;
					t2_8 = 0;	
					for (int m=0; m< nsm_occ; m++){
						double t2_3a = 0;
						double t2_3b = 0;	
						for ( int e=0; e< nsm_vir; e++){
							t2_3a += t_ai(j,e)*gen_Fme(m,e);
							t2_3b += t_ai(i,e)*gen_Fme(m,e);
							//cout << "e from me" << '\t' << e << endl;
						}
						t2_3 += t_aibj(i,m,a,b)*(gen_Fmi(m,j) + 0.5*t2_3a);
						t2_3 -= t_aibj(j,m,a,b)*(gen_Fmi(m,i) + 0.5*t2_3b);
						//cout << "me" << '\t' << m << endl;
						for (int n=0; n < nsm_occ; n++){
							t2_4 += tau(m,n,a,b)*gen_Wmnij(m,n,i,j);
							//cout << "mn" << '\t' << m << '\t' <<  n << endl;
						}
						for (int e=0; e< nsm_vir; e++){
							t2_6 += t_aibj(i,m,a,e)*gen_Wmbej(m,b,e,j) - t_ai(i,e)*t_ai(m,a)*spi_TMO(m,b+nsm_occ,e+nsm_occ,j);
							t2_6 -= t_aibj(i,m,b,e)*gen_Wmbej(m,a,e,j) - t_ai(i,e)*t_ai(m,b)*spi_TMO(m,a+nsm_occ,e+nsm_occ,j);
							t2_6 -= t_aibj(j,m,a,e)*gen_Wmbej(m,b,e,i) - t_ai(j,e)*t_ai(m,a)*spi_TMO(m,b+nsm_occ,e+nsm_occ,i);
							t2_6 += t_aibj(j,m,b,e)*gen_Wmbej(m,a,e,i) - t_ai(j,e)*t_ai(m,b)*spi_TMO(m,a+nsm_occ,e+nsm_occ,i);
							//cout << "me" << '\t' << m<< '\t' << e << endl;
						}
				                t2_8 += t_ai(m,a)*spi_TMO(m,b+nsm_occ,i,j);
						t2_8 -= t_ai(m,b)*spi_TMO(m,a+nsm_occ,a,b);
						//cout << "m" << '\t' << m << endl;
					}
					Dijab = evals(i) + evals(j) - evals(a+nsm_occ) - evals(b+nsm_occ);
                                        t2_1 += spi_TMO(i,j,a+nsm_occ, b+nsm_occ);
					//cout << i << '\t' << j << '\t' << a << '\t' << b << endl;
					update_taibj(i,j,a,b) = t2_1 + t2_2 - t2_3 + t2_4/2 + t2_5/2 + t2_6 + t2_7 - t2_8;
					update_taibj(i,j,a,b) /= Dijab;
					}
	return update_taibj;
}

//Calculate the correlation energy, the first iteration with the initial guess matches with that of MP2 energy

double CCSD::CC_energy(Matrix &t1_amp, Tensor4D &t2_amp){
     	double Ecc = 0;
        double Ecc_1 = 0;
        double Ecc_2 = 0;
        double Ecc_3 = 0;
        for (int i=0; i < nsm_occ; i++)
                for (int a=0; a < nsm_vir; a++){
                        for (int j=0; j < nsm_occ; j++)
                                for (int b=0; b<nsm_vir; b++){
                                        Ecc_2 += spi_TMO(i,j,a+nsm_occ,b+nsm_occ)*t2_amp(i,j,a,b);
                                        Ecc_3 += spi_TMO(i,j,a+nsm_occ,b+nsm_occ)*t1_amp(i,a)*t1_amp(j,b);
                                }
                        Ecc_1 += spi_Focko(i,a+nsm_occ)*t1_amp(i,a);
                }
        Ecc = Ecc_1 + Ecc_2/4 + Ecc_3/2;
	return Ecc;
}

//Convergence condition with energy threshold = 1e-6

bool CCSD::converged(double &energy_dif){
	bool converge = false;
	//commented out the t amplitude threshold 
	/*
	Matrix t1_dif;
	t1_dif.setZero(nsm_occ,nsm_vir);
	Tensor4D t2_dif(nsm_occ, nsm_occ, nsm_vir, nsm_vir);
	t1_dif = new_t1amp - old_t1amp;
	t2_dif = new_t2amp - old_t2amp;
	double t1_norm = t1_dif.norm();
	t2_dif = t2_dif.square();
	Eigen::Tensor<double,0> t2_norm = t2_dif.mean();
	*/
	if (abs(energy_dif) < pow(10,-6)){
		converge = true;
	}
	return converge; 		
}

	

