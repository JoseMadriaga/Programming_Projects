#include "CCSD.h" 

using namespace std;

CCSD CC;

int main() { 

//reads and stores the SCF and MP2 energy

        CC.read_vals("/Users/josemarcmadriaga/Documents/prog_proj/CCSD/Escf.dat", CC.Escf);
        CC.read_vals("/Users/josemarcmadriaga/Documents/prog_proj/CCSD/Emp2.dat", CC.Emp2);
        CC.read_vals("/Users/josemarcmadriaga/Documents/prog_proj/CCSD/Etot.dat", CC.Etot);
        cout << "Escf: ";
        printf("%8.15f\n", CC.Escf);
        cout << "Emp2: ";
        printf("%8.15f\n", CC.Emp2);
        cout << "Etot: ";
        printf("%8.15f\n", CC.Etot);

//static declaration of AO basis 

        CC.nao = 7;
        CC.nocc = 5;
	CC.nvir = CC.nao - CC.nocc;
	CC.nsm = CC.nao*2;
	CC.nsm_occ = CC.nocc*2;
	CC.nsm_vir = CC.nvir*2;

	//cout << CC.nao  << '\t' << CC.nocc << '\t' <<  CC.nvir << '\t' <<  CC.nsm << '\t' <<  CC.nsm_occ << '\t' <<  CC.nsm_vir << endl;

//reads and stores the overlap, kinetic, and nuclear integrals into a matrix
//Same with the AO coefficients

        CC.read_txt("/Users/josemarcmadriaga/Documents/prog_proj/CCSD/s.txt", CC.overlap);
        CC.read_txt("/Users/josemarcmadriaga/Documents/prog_proj/CCSD/t.txt",CC.kinetic);
        CC.read_txt("/Users/josemarcmadriaga/Documents/prog_proj/CCSD/v.txt", CC.nuclear);
        CC.read_Matrix("/Users/josemarcmadriaga/Documents/prog_proj/CCSD/AO_3.dat", CC.AO);

//construct the AO core Hamiltonian

        CC.C_ham.resize(CC.nao, CC.nao);
        CC.C_ham.setZero();
        CC.C_ham = CC.kinetic + CC.nuclear;

//reads and stores the MO two electron integrals in an array 

        CC.read_twoe("/Users/josemarcmadriaga/Documents/prog_proj/CCSD/MO_eric1.dat", CC.MO_eric);

	CC.conv_spi(CC.MO_eric, CC.spi_TMO);

	CC.spi_Fock(CC.AO, CC.C_ham, CC.spi_TMO, CC.spi_Focko);

	CC.guess(CC.spi_Focko, CC.spi_TMO, CC.t_ai, CC.t_aibj);
	
	//CC.tau_and_eff(CC.t_ai, CC.t_aibj, CC.tau_eff, CC.tau);
	
	//CC.up_tai();	

	//CC.up_taibj();
	
//Iteration process of updating the t amplitudes from corresponding intermediates and (effective) taus until the correlation energy converges within a threshold

	double maxiter = 100;
	
	for (int i=1; i < maxiter; i++){
		double old_CCenergy = CC.CC_energy(CC.t_ai,CC.t_aibj);
		CC.tau_and_eff(CC.t_ai, CC.t_aibj, CC.tau_eff, CC.tau);
		CC.up_tai();
		CC.up_taibj();
		double new_CCenergy = CC.CC_energy(CC.update_tai, CC.update_taibj);
		double CCenergy_dif = new_CCenergy - old_CCenergy;
		if (CC.converged(CCenergy_dif)){
			cout << "Converged" << endl;
			printf("Total energy = %6.12f \n", CC.Escf + new_CCenergy);
			break;
		}
		CC.t_ai = CC.update_tai;
		CC.t_aibj = CC.update_taibj;
		cout << "iteration =" << '\t' << i << '\t';
		printf("old_Ecc = %6.12f \t", old_CCenergy);
		printf("new_Ecc = %6.12f \n", new_CCenergy); 
		//cout << CC.t_ai << endl;
	}	

return 0;
}

