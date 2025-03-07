// #include "slip_system/slip_system.h"

// #include <cmath>
// #include <iostream>
// #include <vector>

// namespace slip_system
// {
// double computeShearModulus(double mu_0, double mu_p, double mu_T, double P, double T, double T_r)
// {
//     return mu_0 * (1 + (mu_p / mu_0) * P + (mu_T / mu_0) * (T - T_r));
// }

// double computeDislocationWaitingTime(double v_g, double delta_G, double T)
// {
//     return 1.0 / v_g * (exp(delta_G / (k_B * T)) - 1);
// }

// double computeDislocationRunVelocity(double tau, double tau_ath, double b, double B_T)
// {
//     return (fabs(tau) - tau_ath) * b / B_T * (1.0 - pow((fabs(tau) - tau_ath) * b / (B_T * 0.1), 2));
// }

// double computeDislocationVelocity(double L, double t_w, double t_r)
// {
//     return L / (t_w + t_r);
// }

// double computeShearStrainRate(double rho_m, double b, double v, double tau)
// {
//     return rho_m * b * v * (tau >= 0 ? 1.0 : -1.0);
// }

// double computeMobileDislocationDensity(double rho_m, double rho_im, double gamma_dot, double alpha_het,
//                                        double alpha_mul, double alpha_ann, double alpha_trap, double alpha_depin,
//                                        double dt)
// {
//     double rho_het = alpha_het * fabs(gamma_dot);
//     double rho_mul = alpha_mul * fabs(gamma_dot);
//     double rho_ann = alpha_ann * (rho_m + rho_im) * fabs(gamma_dot);
//     double rho_trap = alpha_trap * sqrt(rho_m) * fabs(gamma_dot);
//     double rho_depin = alpha_depin * sqrt(rho_m) * fabs(gamma_dot);
//     return std::max(rho_m + dt * (rho_het + rho_mul - rho_ann - rho_trap + rho_depin), 0.0);
// }

// double computeSingleShearStrainRate(const std::vector<double>& slip_direction,
//                                     const std::vector<double>& slip_plane_normal, double mu_0, double mu_p, double mu_T,
//                                     double P, double T, double T_r, double b, double v_g, double tau_ath, double rho_im,
//                                     double alpha_het, double alpha_mul, double alpha_ann, double alpha_trap,
//                                     double alpha_depin, double dt)
// {
//     double mu = computeShearModulus(mu_0, mu_p, mu_T, P, T, T_r);

//     // 计算剪切应力 (tau)
//     double tau = 0.0;
//     for (size_t i = 0; i < slip_plane_normal.size(); ++i)
//     {
//         tau += slip_plane_normal[i] * slip_direction[i];
//     }

//     double delta_G = 0.5 * b * b * mu;  // 激活能量 (J)
//     double t_w = computeDislocationWaitingTime(v_g, delta_G, T);

//     double B_T = 1.051e-5 * (1 + 0.9 * (T / 230));  // 温度相关阻尼系数
//     double v_r = computeDislocationRunVelocity(tau, tau_ath, b, B_T);

//     double L = 1.0 / sqrt(rho_im);  // 计算位错间距 (m)
//     double t_r = L / v_r;           // 计算位错运行时间 (s)
//     double v = computeDislocationVelocity(L, t_w, t_r);

//     double gamma_dot = computeShearStrainRate(rho_im, b, v, tau);

//     return gamma_dot;
// }

// std::vector<double> computeShearStrainRates(
//     const std::vector<std::pair<std::vector<double>, std::vector<double>>>& slip_systems, double mu_0, double mu_p,
//     double mu_T, double P, double T, double T_r, double b, double v_g, double tau_ath, double rho_im, double alpha_het,
//     double alpha_mul, double alpha_ann, double alpha_trap, double alpha_depin, double dt)
// {
//     std::vector<double> gamma_dots;

//     for (const auto& system : slip_systems)
//     {
//         const std::vector<double>& slip_direction = system.first;
//         const std::vector<double>& slip_plane_normal = system.second;

//         double gamma_dot =
//             computeSingleShearStrainRate(slip_direction, slip_plane_normal, mu_0, mu_p, mu_T, P, T, T_r, b, v_g,
//                                          tau_ath, rho_im, alpha_het, alpha_mul, alpha_ann, alpha_trap, alpha_depin, dt);

//         gamma_dots.push_back(gamma_dot);
//     }

//     return gamma_dots;
// }
// }  // namespace slip_system