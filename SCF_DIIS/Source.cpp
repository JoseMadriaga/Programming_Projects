#include "SCF.h"
#include "molecule.h"

using namespace std;

SCF fun;

int main() {

	
	Molecule mol("/Users/josemarcmadriaga/Documents/prog_proj/proj_3/geom.dat",0);
	cout << "Numbers of atoms:" << mol.natom << endl;
	cout << "Input Cartersian coordinates:\n";
	mol.print_geom();

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

	// DIIS procedure

	return 0;
}
