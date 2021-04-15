#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using namespace std;

class Molecule
{
public:
    Matrix w;
    Matrix evals;
    double mi;
    double mj;
    Matrix H;
    int natom;
    int charge;
    int* zvals;
    double** geom;
    string point_group;

    void freq();
    void diag_hessian();
    void mass_weight();
    void read_hessian(const char* filename);
    void print_geom();
    void rotate(double phi);
    void translate(double x, double y, double z);
    double bond(int atom1, int atom2);
    double angle(int atom1, int atom2, int atom3);
    double torsion(int atom1, int atom2, int atom3, int atom4);
    double unit(int cart, int atom1, int atom2);
    double oop(int atom1, int atom2, int atom3, int atom4);

    Molecule(const char* filename, int q);
    ~Molecule();
};
