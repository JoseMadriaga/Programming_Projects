#include "CCSD.h"

using namespace std;

double CCSD::gen_taibjck_c(int i, int j, int k, int a, int b, int c){
	double tc = 0;
	for (int e=0; e < nsm_vir; e++){
		tc += t_aibj(j,k,a,e)*spi_TMO(e+nsm_occ, i, b+nsm_occ, c+nsm_occ);  
	}	
	for (int m=0; m < nsm_occ; m++){
		tc -= t_aibj(i,m,b,c)*spi_TMO(m,a+nsm_occ,j,k);
	}
	return tc;
}

double CCSD::gen_taibjck_d(int i, int j, int k, int a, int b, int c){
	double td = 0;
	td = t_ai(i,a)*spi_TMO(j,k,b+nsm_occ,c+nsm_occ);
	return td; 
}

double CCSD::tc_perabc(int i, int j, int k, int a, int b, int c){
	double per_abc = 0;
	per_abc = gen_taibjck_c(i,j,k,a,b,c) - gen_taibjck_c(i,j,k,b,a,c) - gen_taibjck_c(i,j,k,c,b,a);
	return per_abc;
}

double CCSD::tc_perijk(int i, int j, int k, int a, int b, int c){
	double per_ijk = 0;
	per_ijk = tc_perabc(i,j,k,a,b,c) - tc_perabc(j,i,k,a,b,c) - tc_perabc(k,j,i,a,b,c);
	return per_ijk;
}

double CCSD::td_perabc(int i, int j, int k, int a, int b, int c){
        double per_abc = 0;
        per_abc = gen_taibjck_d(i,j,k,a,b,c) - gen_taibjck_d(i,j,k,b,a,c) - gen_taibjck_d(i,j,k,c,b,a);
        return per_abc;
}

double CCSD::td_perijk(int i, int j, int k, int a, int b, int c){
        double per_ijk = 0;
        per_ijk = td_perabc(i,j,k,a,b,c) - td_perabc(j,i,k,a,b,c) - td_perabc(k,j,i,a,b,c);
        return per_ijk;
}

void CCSD::parenT_energy(){
	double parenT_energy = 0;
	double taibjck_c = 0;
	double taibjck_d = 0;
	double Dijkabc = 0;

	for (int i=0; i < nsm_occ; i++)
		for (int j=0; j < nsm_occ; j++)
			for (int k=0; k < nsm_occ; k++)
				for (int a=0; a < nsm_vir; a++)
					for (int b=0; b < nsm_vir; b++)
						for (int c=0; c < nsm_vir; c++){
							taibjck_c = 0;
							taibjck_d = 0;
							Dijkabc = 0;
							taibjck_c = tc_perijk(i,j,k,a,b,c);
							taibjck_d = td_perijk(i,j,k,a,b,c);
							Dijkabc = spi_Focko(i,i) + spi_Focko(j,j) + spi_Focko(k,k) - spi_Focko(a+nsm_occ, a+nsm_occ) -spi_Focko(b+nsm_occ, b+nsm_occ) - spi_Focko(c+nsm_occ, c+nsm_occ);
							taibjck_c /= Dijkabc;
							taibjck_d /= Dijkabc;
							parenT_energy += taibjck_c*Dijkabc*(taibjck_c + taibjck_d);
						}

	parenT_energy /= 36;
	printf("CCSD(T) correction energy = %6.12f \n", parenT_energy);
	printf("Total energy = %6.12f \n", CCSD_energy + parenT_energy);
} 
