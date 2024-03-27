// Port of 2F1 code from https://arxiv.org/abs/0708.0116v2
// Authors: N. Michel, M. V. Stoitsov
//
// ---------------------------------------------------------------------------

#ifndef HYP_2F1
#define HYP_2F1

#include "constants.hpp"

namespace analyticRT
{
    inline double inf_norm (const complex &z)
    {
        return std::max(std::abs(std::real (z)), std::abs(std::imag (z)));
    };

    // Test of finiteness of a complex number
    // --------------------------------------
    // If real or imaginary parts are finite, true is returned.
    // Otherwise, false is returned

    inline bool isfinite (const complex &z)
    {
        const double x = std::real(z), y = std::imag(z);
        return (finite(x) && finite(y));
    }

    complex log1p (const complex &z);
    complex expm1 (const complex &z);

    complex Gamma_inv (const complex &z);
    complex log_Gamma (const complex &z);

    complex Gamma_ratio_diff_small_eps (const complex &z,const complex &eps);
    complex Gamma_inv_diff_eps (const complex &z,const complex &eps);
    complex A_sum_init (const int m,const complex &eps,const complex &Gamma_inv_one_meps);
    complex log_A_sum_init (const int m,const complex &eps);
    complex B_sum_init_PS_one (const complex &a,const complex &b,const complex &c,
                    const complex &Gamma_c,const complex &Gamma_inv_one_meps,
                    const complex &Gamma_inv_eps_pa_pm,const complex &Gamma_inv_eps_pb_pm,
                    const complex &one_minus_z,const int m,const complex &eps);
    complex B_sum_init_PS_infinity (const complex &a,const complex &c,
                        const complex &Gamma_c,const complex &Gamma_inv_cma,
                        const complex &Gamma_inv_one_meps,const complex &Gamma_inv_eps_pa_pm,
                        const complex &z,const int m,const complex &eps);
    void cv_poly_der_tab_calc (const complex &a,const complex &b,const complex &c,const complex &z,double cv_poly_der_tab[]);
    double cv_poly_der_calc (const double cv_poly_der_tab[],const double x);
    int min_n_calc (const double cv_poly_der_tab[]);
    complex hyp_PS_zero (const complex &a,const complex &b,const complex &c,const complex &z);
    complex hyp_PS_one (const complex &a,const complex &b,const complex &c,const complex &one_minus_z);
    complex hyp_PS_infinity (const complex &a,const complex &b,const complex &c,const complex &z);
    complex hyp_PS_complex_plane_rest (const complex &a,const complex &b,const complex &c,const complex &z);
    complex hyp_2F1 (const complex &a,const complex &b,const complex &c,const complex &z);
};
#endif