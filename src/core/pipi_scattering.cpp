// GKPY [1] parameterizations for pi pi phase-shifts and inelasticities
// 
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------
// REFERENCES 
// [1] - https://arxiv.org/abs/1102.2183
// ---------------------------------------------------------------------------

#include "pipi_scattering.hpp"

namespace analyticRT
{
    //-----------------------------------------------------------------------------
    // Produces the Phase-shift for given Partial Wave as function of s.
    // Based on parameterizations from CFD Set:
    // R. Garcia-Martin, R. Kaminski, J. R. Palaez, J. Ruiz de Elvira, F. J. Yndurain
    // Phys. Rev. D 83, 074004 (2011).
    //
    // Ported from A. Jackura and N. Sherrill's fortran code
    // ajackura@indiana.edu, nlsherri@indiana.edu
    double pipi::phase_shift(int iso, int l, double s)
    {
        //Error Checks :)
        int wave = (3*l - iso)/2;
        if (iso > 2 || iso < 0) {wave = 600;}
        if (l % 2 != iso % 2) {wave = 600;}
        if (l >= 4 && wave != 600)
        {std::cout << "Phase shift contributions of j > 3 are completely negligible. Quitting... \n";
        exit(1);}

        double cot_delta, delta, sh;
        double temp1, temp2, temp3;
        sh = pow(1.42, 2.);

        //Momenta
        double k, k2;
        k = elastic_mom(s, 4.*M_PION*M_PION);
        k2 = elastic_mom(s, 4.*M_KAON*M_KAON);

        //S2 wave (l = 0, iso = 2)
        if (wave == -1)
        {
            double B0, B1, z2, wl;
            double sl, sm;

            sl = pow(1.05, 2.);
            sm = pow(.850, 2.);

            wl = conformal(s, sl);

            z2 = .1435;
            B0 = -79.4;
            B1 = -63.0;

            temp1 = sqrt(s) / (2.*k);
            temp2 = M_PION*M_PION / (s - 2.*z2*z2);

            if ((s > 4.*M_PION*M_PION) && (s <= sm)) //Low energy parameterization
            {
                temp3 = B0 + B1 * wl;
                cot_delta = temp1 * temp2 * temp3;
                delta = atan2(1., cot_delta);
            }
            else if ((s > sm) && ( s < sh)) //Intermediate energies
            {
                double Bh0, Bh1, Bh2, wh, whsm, wlsm;
                double temp4;

                wh = conformal(s, sh);
                whsm = conformal(sm, sh);
                wlsm = conformal(sm, sl);

                Bh2 = 32.;
                Bh0 = B0 + B1 * wlsm;
                temp4 = pow((sqrt(sm) + sqrt(sh - sm)) / (sqrt(sm) + sqrt(sl - sm)), 2.);
                Bh1 = B1 * (sl / sh) * (sqrt(sh - sm) / sqrt(sl - sm)) * temp4;

                temp3 = Bh0 + Bh1*(wh - whsm) + Bh2*pow(wh - whsm, 2.);

                cot_delta = temp1 * temp2 * temp3;
                delta = atan2(1., cot_delta);
            }
            else {delta = 0.;}

            return delta;
        }

        //S0 wave (l = 0, iso = 0)
        if (wave == 0)
        {
            double sm = pow(0.85, 2.);
            double B0, B1, B2, B3, w;
            B0 = 7.14;
            B1 = -25.3;
            B2 = -33.2;
            B3 = -26.2;
            w = conformal(s, 4.*M_KAON*M_KAON);

            if ((s > 4.*M_PION*M_PION) && (s <= sm)) //Low energy parameterization
            {
                temp1 = sqrt(s) / (2.*k);
                temp2 = M_PION * M_PION / (s - 0.5 * M_PION * M_PION);
                temp3 = M_PION / sqrt(s);

                cot_delta = temp1 * temp2 * (temp3 + B0 + B1 * w + B2 * w * w + B3 *w*w*w);
                delta = atan2(1., cot_delta);
            }
            else if (s < sh)
            {
                double d0, C1, B, C2, D;
                double k2, k2m;

                d0 = 226.5 * DEG2RAD; //DEG2RADerting to radians
                C1 = -81. * DEG2RAD;
                B = 93.3 * DEG2RAD;
                C2 = 48.7 * DEG2RAD;
                D = -88.3 * DEG2RAD;

                k2m = elastic_mom(4.*M_KAON*M_KAON, sm); //Switched inputs to avoid negative under radical

                if ((s > sm) && (s < 4.*M_KAON*M_KAON)) // intermediate energy
                {
                    double cot_delm, delm, delPm, km, wm;
                    double temp4, temp5, temp6;
                    k2 = elastic_mom(4.*M_KAON*M_KAON, s);
                    temp4 = pow(1. - (k2 / k2m), 2.);
                    temp5 = k2 * (2. - (k2 / k2m)) / k2m;
                    temp6 = (k2m - k2) / pow(M_KAON, 3.);


                    //Calculate delta(sm)
                    wm = conformal(sm, 4.*M_KAON*M_KAON);
                    km = elastic_mom(sm, 4.*M_PION*M_PION);
                    temp1 = sqrt(sm) / (2. * km);
                    temp2 = M_PION * M_PION / (sm - 0.5 * M_PION * M_PION);
                    temp3 = M_PION / sqrt(sm);
                    cot_delm = temp1 * temp2 * (temp3 + B0 + B1 * wm + B2 * wm * wm + B3 *wm*wm*wm);
                    delm = atan2(1., cot_delm);

                    //Derivative calculated in Mathematica (because im lazy)
                    delPm = 1.588503;

                    delta = d0 * temp4 + delm * temp5 + k2 * (k2m - k2) * (8.*delPm + C1 * temp6);
                }

                if ((s > 4.*M_KAON*M_KAON) && (s < sh)) //above KK threshold
                {
                    k2 = elastic_mom(s, 4.*M_KAON*M_KAON);
                    temp1 = (k2*k2) / (M_KAON * M_KAON);
                    delta = d0 + B*temp1 + C2 * temp1*temp1;

                    if (s > 4.*M_ETA*M_ETA)
                    {
                        double k3 = elastic_mom(s, 4.*M_ETA*M_ETA);
                        temp2 = D * (k3*k3) / (M_ETA * M_ETA);
                        delta += temp2;
                    }
                }
            }
            return delta;
        }

    //P1 wave (l = 1, iso = 1)
        else if (wave == 1)
        {
            double s0, w;

            s0 = pow(1.05, 2.);
            w = conformal(s, s0);

            if ((s > 4.*M_PION*M_PION) && (s < 4.*M_KAON*M_KAON))
            {
                double B0, B1;
                B0 = 1.043;
                B1 = .19;

                temp1 = sqrt(s)*(M_RHO*M_RHO - s)/(2. * pow(k, 3.));
                temp2 = 2.*pow(M_PION, 3.)/(M_RHO*M_RHO*sqrt(s));

                cot_delta = temp1*(temp2 + B0 + B1* w);
                delta = atan2(1., cot_delta);
            }
            else if ((s >= 4.*M_KAON*M_KAON) && (s < sh))
            {
                double lambda0, wK, KK;
                double lambda1, lambda2;
                lambda1 = 1.38;
                lambda2 = -1.70;

                //lambda0 = low energy (see above) at KK threshold
                double B0, B1;
                B0 = 1.043;
                B1 = .19;
                KK = elastic_mom(4.*M_KAON*M_KAON, 4.*M_PION*M_PION);
                wK = conformal(4.*M_KAON*M_KAON, s0);
                temp1 = sqrt(4.*M_KAON*M_KAON)*(M_RHO*M_RHO - s)/(2. * pow(KK, 3.));
                temp2 = 2.*pow(M_PION, 3.)/(M_RHO*M_RHO*sqrt(4.*M_KAON*M_KAON));
                lambda0 = temp1*(temp2 + B0 + B1* wK);

                temp3 = (sqrt(s) / (2.*M_KAON)) - 1.;

                cot_delta = lambda0 + lambda1*temp3 + lambda2*temp3*temp3;
                delta = atan2(1., cot_delta);
            }
            else {delta = 0.;}
            return delta;
        }

        //D2 wave (l = 2, iso = 2)
        else if (wave == 2)
        {
            if ((s > 4.*M_PION*M_PION) && (s < sh))
            {
                double w, s0, Del;
                double B0, B1, B2;

                s0 = pow(1.45,2.);
                B0 = 4.1e3;
                B1 = 8.6e3;
                B2 = 25.5e3;
                Del = .233;

                w = conformal(s, s0);

                temp1 = sqrt(2) / (2.* pow(k, 5.));
                temp2 = (pow(M_PION, 4.) * s) / (4.*M_PION*M_PION + 4.*Del*Del -s);
                temp3 = B0 + B1 * w + B2 * w * w;

                cot_delta = temp1 * temp2 * temp3;
                delta = atan2(1., cot_delta);
            }
            else {delta = 0.;}
            return delta;
        }

        //D0 wave (l = 2, iso = 0)
        else if (wave == 3)
        {
            double B0, B1;
            double s0 = pow(1.05, 2.);

            B0 = 12.40;
            B1 = 10.06;

            temp1 = sqrt(s) / (2.*pow(k, 5.));
            temp2 = (M_F2*M_F2 - s)*M_PION*M_PION;

            if ((s > 4.*M_PION*M_PION) && (s <= 4.*M_KAON*M_KAON))
            {
                double w = conformal(s, s0);
                temp3 = B0 + B1* w;
                cot_delta = temp1 * temp2 * temp3;
                delta = atan2(1., cot_delta);
            }
            else if ((s > 4.*M_KAON*M_KAON) && (s < sh))
            {
                double wh, s02, B0h, B1h;

                B1h = 43.2;
                s02 = pow(1.45, 2.);
                wh = conformal(s, s02);
                B0h = 18.69;

                temp3 = B0h + B1h * wh;
                cot_delta = temp1 * temp2 * temp3;
                delta = atan2(1., cot_delta);
            }
            else {delta = 0.;}
            return delta;
        }

        //F1 wave (l = 3, iso = 1)
        else if (wave == 4)
        {
            if ((s > 4.*M_PION*M_PION) && (s < sh))
            {
                double B0, B1, lambda;
                double s0, w;

                s0 = pow(1.45, 2.);
                w = conformal(s, s0);
                B0 = 1.09e5;
                B1 = 1.41e5;
                lambda = 0.051e5;

                temp1 = sqrt(s)*pow(M_PION, 6.) / (2.*pow(k, 7.));
                temp2 = 2.*lambda * M_PION / sqrt(s);

                cot_delta = temp1 * (temp2 + B0 + B1 * w);
                return atan2(1., cot_delta);
            }
        }

        //Unphysical partial waves (ie anything else)
        else
        {
            std::cout << "Invalid Phase Shift with l = " << l << " and isospin = " << iso << ". Quiting... \n";
            exit(1);
        }
        return 0.;
    };

    //-----------------------------------------------------------------------------
    // Produces the Inelasticity for given Partial Wave as a function of s.
    // Based on CFD parameterizations from:
    // R. Garcia-Martin, R. Kaminski, J. R. Palaez, J. Ruiz de Elvira, F. J. Yndurain
    // Phys. Rev. D 83, 074004 (2011).
    //
    // Ported from A. Jackura and N. Sherrill's fortran code
    // ajackura@indiana.edu, nlsherri@indiana.edu
    double pipi::inelasticity(int iso, int l, double s)
    {
        double eta;
        double sh = pow(1.42, 2.);

        int wave = (3*l - iso)/2;
        if (iso > 2 || iso < 0 || l < 0) {wave = 600;}
        if (l % 2 != iso % 2) {wave = 600;} //Check Bose symmetry
        if (l >= 3 && wave != 600) {return eta = 1.;} //Inelasticities for high partial waves are negligable

//S2 wave (l = 0, iso = 2)
        if (wave == -1)
        {
            double si = pow(1.05, 2.);

            if ((s > 4.*M_PION*M_PION) && (s <= si))          //Low energy parameterization
            {
                eta = 1;
            }
            else if ((s > si) && ( s <= sh))        //Intermediate energies
            {
                double eps = .28;

                eta = 1. - eps*pow(1 - (si/s), 1.5);
            }
            return eta;
        }

    //S0 wave (l = 0, iso = 0)
        else if (wave == 0)
        {
            if ((s > 4.*M_PION*M_PION) && (s <= 4.*M_KAON*M_KAON))
            {
                    eta = 1.;
            }
            else if ((s > 4.*M_KAON*M_KAON) && (s < sh))
            {
                double ep1, ep2, ep3;
                double k2;
                double temp1, temp2, temp3;

                ep1 = 4.9;
                ep2 = -15.1;
                ep3 = 4.7;
                k2 = elastic_mom(s, 4.*M_KAON*M_KAON);

                temp1 = -k2 / sqrt(s);
                temp2 = ep1 + ep2 * (k2/sqrt(s)) + ep3 * (k2 * k2 / s);
                temp3 = temp1 * temp2 * temp2;
                if (s > 4.*M_ETA*M_ETA)
                {
                    double ep4, k3;
                    ep4 = 0.32;
                    k3 = elastic_mom(s, 4.*M_ETA*M_ETA);

                    temp3 += -ep4 * k3 / sqrt(s);
                }
                eta = exp(temp3);
            }
            return eta;
        }

    //P1 wave (l = 1, iso = 1)
        else if (wave == 1)
        {
            if ((s > 4.*M_PION*M_PION) && (s < 4.*M_KAON*M_KAON))
            {
                eta = 1.;
            }
            else if ((s > 4.*M_KAON*M_KAON) && (s < sh))
            {
                double temp;
                double ep1, ep2;
                ep1 = 0.;
                ep2 = .07;

                temp = 1. - (4.*M_KAON*M_KAON/s);
                eta = 1 - ep1 * sqrt(temp) - ep2 * temp;
            }
            return eta;
        }

    //D2 wave (l = 2, iso = 2)
        else if (wave == 2)
        {
            double si = pow(1.05, 2.);
            if ((s > 4.*M_PION*M_PION) && (s <= si))
            {
                eta = 1.;
            }
            else if ((s > si) && (s < sh))
            {
                double ep = 0.;
                eta = 1. - ep * pow(1 - (si/s),3.);
            }
            return eta;
        }

    //D0 wave (l = 2, iso = 0)
        else if (wave == 3)
        {
            if ((s > 4.*M_PION*M_PION) && (s <= 4.*M_KAON*M_KAON))
            {
                eta = 1.;
            }
            if ((s > 4.*M_KAON*M_KAON) && (s < sh))
            {
                double eps, r, temp1, temp2;
                eps = 0.254;
                r = 2.29;

                double k2 = elastic_mom(s, 4.*M_KAON*M_KAON);

                temp1 = (1. - 4.*M_KAON*M_KAON/s) / (1. - 4.*M_KAON*M_KAON/(M_F2*M_F2));
                temp2 = 1. - k2 / elastic_mom(M_F2*M_F2, 4.*M_KAON*M_KAON);

                eta = 1. - eps*pow(temp1, 2.5)*(1. + r*temp2);
            }
            return eta;
        }

    //Unphysical partial waves (ie anything else)
        else
        {
            std::cout << "Invalid Inelasticity with l = " << l << " and isospin = " << iso << ". Quiting... \n";
            exit(1);
        }
    }
};