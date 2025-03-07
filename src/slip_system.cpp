#include "slip_system/slip_system.h"

namespace slip_system
{
// 计算剪切模量
double computeShearModulus(double mu_0, double mu_p, double mu_T, double P, double T, double T_r)
{
    return mu_0 * (1 + (mu_p / mu_0) * P + (mu_T / mu_0) * (T - T_r));
}

// 计算位错等待时间
double computeDislocationWaitingTime(double v_g, double delta_G, double T)
{
    return 1.0 / v_g * (exp(delta_G / (k_B * T)) - 1);
}

// 计算位错运行速度
double computeDislocationRunVelocity(double tau, double tau_ath, double b, double B_T)
{
    return (fabs(tau) - tau_ath) * b / B_T * (1.0 - pow((fabs(tau) - tau_ath) * b / (B_T * 0.1), 2));
}

// 计算位错速度
double computeDislocationVelocity(double L, double t_w, double t_r)
{
    return L / (t_w + t_r);
}

// 计算剪切应变率
double computeShearStrainRate(double rho_m, double b, double v, double tau)
{
    return rho_m * b * v * (tau >= 0 ? 1.0 : -1.0);
}

// 计算可动位错密度
double computeMobileDislocationDensity(double rho_m, double rho_im, double gamma_dot, double alpha_het,
                                       double alpha_mul, double alpha_ann, double alpha_trap, double alpha_depin,
                                       double dt)
{
    double rho_het = alpha_het * fabs(gamma_dot);
    double rho_mul = alpha_mul * fabs(gamma_dot);
    double rho_ann = alpha_ann * (rho_m + rho_im) * fabs(gamma_dot);
    double rho_trap = alpha_trap * sqrt(rho_m) * fabs(gamma_dot);
    double rho_depin = alpha_depin * sqrt(rho_m) * fabs(gamma_dot);
    return std::max(rho_m + dt * (rho_het + rho_mul - rho_ann - rho_trap + rho_depin), 0.0);
}

}  // namespace slip_system