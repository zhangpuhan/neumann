//
// Created by Puhan Zhang on 1/12/21.
//

#include "model.h"
#include "util.h"


void Square_Lattice::init_lattice() {

    for (int x = 0; x < L; x++)
        for (int y = 0; y < L; y++) {
            int idx = index(x, y);
            site[idx].idx = idx;
            site[idx].x = x;
            site[idx].y = y;

            site[idx].sgn = ((x + y) % 2 == 0) ? +1 : -1;

        }

    for (int i = 0; i < Ns; i++) {
        int j;
        int x = site[i].x;
        int y = site[i].y;

        j = index(mod(x + 1, L), y);
        site[i].nn1[0] = &site[j];

        j = index(mod(x - 1, L), y);
        site[i].nn1[1] = &site[j];

        j = index(x, mod(y + 1, L));
        site[i].nn1[2] = &site[j];


        j = index(x, mod(y - 1, L));
        site[i].nn1[3] = &site[j];

    }
}

arma::sp_cx_mat Square_Lattice::build_Hamiltonian() {

    arma::sp_cx_mat H(dim, dim);
    for(int i = 0; i < Ns; i++) {
        for(auto & k : site[i].nn1) {
            int j = k->idx;
            H(i, j) = t1;
        }
    }

    for(int i = 0; i < Ns; i++) {
        H(i, i) += onsite_V[i];
    }

    return H;
}

void Square_Lattice::integrate_EOM_RK4(double dt) {
    
    arma::sp_cx_mat H(Hamiltonian);
    arma::cx_mat D(Density_Mat);
    arma::cx_mat D2, KD, KD_sum;

    
    // ------- RK4 step-1: ----------------
    
    KD = -_I * dt * ( H * D - D * H );
    
    D2 = D + 0.5 * KD;
    KD_sum = KD / 6.;
    
    
    // ------- RK4 step-2: ----------------
    
    H = build_Hamiltonian();
    
    KD = -_I * dt * ( H * D2 - D2 * H );

    
    D2 = D + 0.5 * KD;
    KD_sum += KD / 3.;

    
    // ------- RK4 step-3: ----------------
    
    H = build_Hamiltonian();
    KD = -_I * dt * ( H * D2 - D2 * H );
    
    
    D2 = D + KD;
    KD_sum += KD / 3.;

    
    // ------- RK4 step-4: ----------------
    
    H = build_Hamiltonian();
    KD = -_I * dt * ( H * D2 - D2 * H );
    KD_sum += KD / 6.;
    
    
    // ------- RK4: sum all steps: ------------
    
    Density_Mat = D + KD_sum;
    
    // compute the system Hamiltonian, R, Delta:
    
    Hamiltonian = build_Hamiltonian();
    
}

void Square_Lattice::compute_fermi_level(arma::vec &eigE) {

    double x1 = eigE(0);
    double x2 = eigE(eigE.size() - 1);

    int max_bisection = 500;
    double eps_bisection = 1.e-12;

    int iter = 0;
    while(iter < max_bisection || fabs(x2 - x1) > eps_bisection) {

        double xm = 0.5 * (x1 + x2);
        double density = 0;
        for(int i=0; i<eigE.size(); i++) {
            density += fermi_density(eigE(i), kT, xm);
        }
        density /= ((double) dim);

        if(density <= filling) x1 = xm;
        else x2 = xm;

        iter++;
    }

    mu = 0.5*(x1 + x2);
}

arma::cx_mat Square_Lattice::compute_density_matrix(arma::sp_cx_mat & H) {

    arma::cx_mat Hm(H);
    arma::cx_mat D(dim, dim);

    arma::vec eigval(dim);
    arma::cx_mat eigvec;

    arma::eig_sym(eigval, eigvec, Hm);

    std::ofstream fs("eig.dat");
    for(int r=0; r<dim; r++) {
        fs << eigval(r) << std::endl;
    }
    fs.close();

    compute_fermi_level(eigval);

    arma::vec fd_factor(dim);
    for(int i=0; i<dim; i++) fd_factor(i) = fermi_density(eigval(i), kT, mu);

    D.zeros();
    for(int a=0; a<dim; a++)
        for(int b=a; b<dim; b++) {

            arma::cx_double sum = 0;
            for(int m=0; m<dim; m++) {
                sum += fd_factor(m) * conj(eigvec(a, m)) * eigvec(b, m);
            }
            D(a, b) = sum;
            if(a != b) D(b, a) = conj(sum);
        }

    return D;
}

void Square_Lattice::save_configuration(std::string const filename) {

    std::ofstream fs;

    fs.open(filename.c_str(), std::ios::out);
    fs.precision(12);

    for(int i=0; i<Ns; i++) {

        fs << real(Density_Mat(2 * i, 2 * i) + Density_Mat(2 * i + 1, 2 * i + 1)) << '\t';
        fs << std::endl;
    }
    fs.close();
}

void Square_Lattice::init_quenched_disorder(double W) {

    std::random_device seed;

    RNG rng = RNG(seed());

    std::uniform_real_distribution<double> rd(-W, W);

    for(int i=0; i<Ns; i++) onsite_V[i] = rd(rng);
}

void Square_Lattice::simulate_dynamics(int max_steps, double dt, double W) {
    // init_quenched_disorder(W);
    Hamiltonian = build_Hamiltonian();
    Density_Mat = compute_density_matrix(Hamiltonian);
//    std::cout << Hamiltonian << std::endl;
//    std::cout << Density_Mat << std::endl;

    time = 0;
    for(int i = 0; i < max_steps; i++) {

        std::cout << "i = " << i << std::endl;

        integrate_EOM_RK4(dt);
        time += dt;
        if (time >= 0.1) {
            init_quenched_disorder(0.1);
        }
        double n_tot = real( arma::trace(Density_Mat));
        std::cout << time << "," << n_tot << std::endl;

    }

}
