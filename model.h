//
// Created by Puhan Zhang on 1/12/21.
//

#ifndef NEUMANN_MODEL_H
#define NEUMANN_MODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <random>
#include "armadillo"

typedef std::mt19937 RNG;

const arma::cx_double _I(0, 1);

class Square_Lattice {
public:
    int L, Ns;
    int dim;

    double t1;

    double filling;
    double mu;
    double kT;

    static constexpr int N_nn1 = 4;

    double time;

    class Site {
    public:
        int idx;
        int x, y;

        int sgn;

        Site *nn1[N_nn1];

    } *site;

    double *onsite_V;

    arma::sp_cx_mat Hamiltonian;
    arma::cx_mat Density_Mat;

    Square_Lattice(int linear_size) {

        L = linear_size;
        Ns = L * L;
        dim = Ns;

        onsite_V = new double[Ns];

        site = new Site[Ns];

        Hamiltonian = arma::sp_cx_mat(dim, dim);
        Density_Mat = arma::cx_mat(dim, dim);

        init_lattice();

        time = 0;

    };

    ~Square_Lattice() {

        delete [] site;
        site = nullptr;
    };

    inline int index(int x, int y) {
        return L * y + x;
    };

    void init_lattice();

    arma::sp_cx_mat build_Hamiltonian();

    void compute_fermi_level(arma::vec &eigE);
    arma::cx_mat compute_density_matrix(arma::sp_cx_mat &);

    // ========================================

    void save_configuration(std::string const filename);
    void simulate_dynamics(int max_steps, double dt, double W);
    void init_quenched_disorder(double W);
    void integrate_EOM_RK4(double dt);



};

#endif //NEUMANN_MODEL_H
