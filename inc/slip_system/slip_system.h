#ifndef SLIP_SYSTEM_H
#define SLIP_SYSTEM_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace slip_system
{

const int N_S = 12;           // 滑移系统数量
const double k_B = 1.38e-23;  // 玻尔兹曼常数 (J/K)

// 计算剪切模量
double computeShearModulus(double mu_0, double mu_p, double mu_T, double P, double T, double T_r);

// 计算位错等待时间
double computeDislocationWaitingTime(double v_g, double delta_G, double T);

// 计算位错运行速度
double computeDislocationRunVelocity(double tau, double tau_ath, double b, double B_T);

// 计算位错速度
double computeDislocationVelocity(double L, double t_w, double t_r);
// 计算剪切应变率
double computeShearStrainRate(double rho_m, double b, double v, double tau);

// 计算可动位错密度
double computeMobileDislocationDensity(double rho_m, double rho_im, double gamma_dot, double alpha_het,
                                       double alpha_mul, double alpha_ann, double alpha_trap, double alpha_depin,
                                       double dt);

}  // namespace slip_system

#endif  // !SLIP_SYSTEM_H
