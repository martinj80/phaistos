// fb5.cpp --- Density and sampler for the FB5 distribution
// Copyright (C) 2006-2008 Wouter Boomsma, Thomas Hamelryck
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

#include "fb5.h"

namespace phaistos {

/////////////
// Density //
/////////////

// Unit vector to array of theta/tau
void Fb5Density::xyz_to_polar(Vector_3D *x, double *tt) {

    tt[0]=acos((*x)[0]);
    tt[1]=atan2((*x)[2], (*x)[1]);
}

// Unit theta/tau angles to x,y,z vector
void Fb5Density::polar_to_xyz(double theta, double tau, Vector_3D *v) {
    double x[3];

    x[0]=cos(theta);
    x[1]=sin(theta)*cos(tau);
    x[2]=sin(theta)*sin(tau);

    (*v)=x;
}

// Constructor
Fb5Density::Fb5Density(double kappa, double beta,
                       std::vector<double> &e1, std::vector<double> &e2, std::vector<double> &e3) {
     this->kappa = kappa;
     this->beta = beta;;
     this->e1 = new Vector_3D(e1);
     this->e2 = new Vector_3D(e2);
     this->e3 = new Vector_3D(e3);

    // Normalizing constant
    logc=std::log(2*M_PI)+kappa-std::log(sqrt((kappa-2*beta)*(kappa+2*beta)));
}

// Copy constructor
Fb5Density::Fb5Density(const Fb5Density &density) {
     kappa = density.kappa;
     beta = density.beta;
     logc = density.logc;
     e1 = new Vector_3D(*density.e1);
     e2 = new Vector_3D(*density.e2);
     e3 = new Vector_3D(*density.e3);
}

// Destructor
Fb5Density::~Fb5Density() {
    delete e1;
    delete e2;
    delete e3;
}




/////////////
// Sampler //
/////////////


// Constructor
Fb5Sampler::Fb5Sampler(double kappa, double beta,
                       std::vector<double> &e1, std::vector<double> &e2, std::vector<double> &e3,
                       RandomNumberEngine *random_number_engine)
     : variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

    a=4*(kappa-2*beta);
    b=4*(kappa+2*beta);
    gamma=(b-a)/2;

    if(b>0) {
        c2=0.5*b/(b-gamma);
    } else {
        c2=0.0;
    }

    lam1=sqrt(a+2*sqrt(gamma));
    lam2=sqrt(b);

    // Put e1,e2,e3 in a 3x3 matrix
    Matrix_3D rho_t(e1, e2, e3);
    rho=new Matrix_3D();

    *rho=transpose(rho_t);
}

// Copy constructor
Fb5Sampler::Fb5Sampler(const Fb5Sampler &other):
     variate_generator_uniform(other.variate_generator_uniform) {

     a = other.a;
     b = other.b;
     c2 = other.c2;
     gamma = other.gamma;
     lam1 = other.lam1;
     lam2 = other.lam2;
     rho = new Matrix_3D(*other.rho);
}

// Copy constructor
Fb5Sampler::Fb5Sampler(const Fb5Sampler &sampler,
                       RandomNumberEngine *random_number_engine)
     : variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

     a = sampler.a;
     b = sampler.b;
     c2 = sampler.c2;
     gamma = sampler.gamma;
     lam1 = sampler.lam1;
     lam2 = sampler.lam2;
     rho = new Matrix_3D(*sampler.rho);
}

// Return next sample from Kent distribution with given parameters
Vector_nD Fb5Sampler::sample() {

    double tt[2];
     
    while(1) {
        double v1, v2, u1, u2;
        double x1, x2, x, y, z;
        double ratio1, ratio2, sign1, sign2;
        double ct, st, cp, sp, phi, theta;

        v1=variate_generator_uniform();
        v2=variate_generator_uniform();
        u1=variate_generator_uniform();
        u2=variate_generator_uniform();

        try {
            x1=-std::log(1-v1*(1-exp(-lam1)))/lam1;
        } catch(...) {
            x1=v1;
        }

        try {
            x2=-std::log(1-v2*(1-exp(-lam2)))/lam2;
        } catch(...) {
            x2=v2;
        }

        if ((x1*x1+x2*x2)>1) {
            continue;
        }

        ratio1=exp(-0.5*(a*x1*x1+gamma*x1*x1*x1*x1)-1+lam1*x1);

        if (u1>ratio1) {
            continue;
        }

        ratio2=exp(-0.5*(b*x2*x2-gamma*x2*x2*x2*x2)-c2+lam2*x2);

        if (u2>ratio2) {
            continue;
        }

        sign1= (variate_generator_uniform()-0.5) < 0 ? -1.0 : 1.0;
        sign2= (variate_generator_uniform()-0.5) < 0 ? -1.0 : 1.0;

        x1=x1*sign1;
        x2=x2*sign2;
        theta=acos(1-2*(x1*x1+x2*x2));
        phi=atan2(x2, x1);
        ct=cos(theta);
        st=sin(theta);
        cp=cos(phi);
        sp=sin(phi);
        x=ct;
        y=st*cp;
        z=st*sp;

        // Return sample as vector
        Vector_3D p(x, y, z);

        // Final sample is matrixmultiply(rho,p)
        // which rotates sample to mean direction
        p=(*rho)*p;

        Fb5Density::xyz_to_polar(&p, tt);

        break;
    }
    return Vector_nD(2, tt);
}


/////////////////////
// Fb5Distribution //
/////////////////////

// Constructor - Initialize sampler and density
Fb5Distribution::Fb5Distribution(double kappa, double beta,
                                 std::vector<double> &e1, 
                                 std::vector<double> &e2, 
                                 std::vector<double> &e3,
                                 RandomNumberEngine *random_number_engine) {

     this->density = new Fb5Density(kappa, beta, e1, e2, e3);
     this->sampler = new Fb5Sampler(kappa, beta, e1, e2, e3, random_number_engine);
}

// Copy constructor
Fb5Distribution::Fb5Distribution(const Fb5Distribution &other) {

     this->density = new Fb5Density(*other.density);
     this->sampler = new Fb5Sampler(*other.sampler);
}

// Copy constructor
Fb5Distribution::Fb5Distribution(const Fb5Distribution &other,
                                 RandomNumberEngine *random_number_engine) {

     this->density = new Fb5Density(*other.density);
     this->sampler = new Fb5Sampler(*other.sampler, random_number_engine);
}

// Destructor
Fb5Distribution::~Fb5Distribution() {
     delete density;
     delete sampler;
}

// Draw sample
Vector_nD Fb5Distribution::sample() {
     return sampler->sample();
}

// Evaluate log-likelihood
double Fb5Distribution::get_log_likelihood(double *angle_pair) {
     return density->get_log_likelihood(angle_pair);
}


}
