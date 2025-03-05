#ifndef FORMULA_H

#include <cmath>  // 用于数学计算
#include <iostream>
#include <vector>

namespace formula
{

// 物理常数
const double k_B = 1.38e-23;  // 玻尔兹曼常数 (J/K)
const double C_t = 5000.0;    // 剪切波速 (m/s)

// 计算剪切模量 (受温度影响)
// mu_0: 初始剪切模量 (Pa)
// mu_p: 剪切模量对压力的影响系数
// mu_T: 剪切模量对温度的影响系数
// P: 压力 (Pa)
// T: 当前温度 (K)
// T_r: 参考温度 (K)
double computeShearModulus(double mu_0, double mu_p, double mu_T, double P, double T, double T_r);

// 计算位错的等待时间 (热激活效应控制)
// v_g: 位错的尝试频率 (s^-1)
// delta_G: 激活能量 (J)
// T: 温度 (K)
double computeDislocationWaitingTime(double v_g, double delta_G, double T);

// 计算位错的运行速度 (受剪应力驱动)
// tau: 剪应力 (Pa)
// tau_ath: 无热临界应力 (Pa)
// b: Burgers 矢量 (m)
// B_T: 温度相关的阻尼系数 (Pa·s)
double computeDislocationRunVelocity(double tau, double tau_ath, double b, double B_T);

// 计算位错速度 (结合等待时间和运行时间)
// L: 位错间距 (m)
// t_w: 位错等待时间 (s)
// t_r: 位错运行时间 (s)
double computeDislocationVelocity(double L, double t_w, double t_r);

// 计算剪切应变率 (Orowan 公式)
// rho_m: 可动位错密度 (m^-2)
// b: Burgers 矢量 (m)
// v: 位错平均速度 (m/s)
// tau: 剪应力 (Pa)
double computeShearStrainRate(double rho_m, double b, double v, double tau);

// 计算可动位错密度 (rho_m)，基于位错的成核、增殖、湮灭等机制
double computeMobileDislocationDensity(double rho_m, double rho_im, double gamma_dot, double alpha_het,
                                       double alpha_mul, double alpha_ann, double alpha_trap, double alpha_depin,
                                       double dt);

}  // namespace formula

#endif  // !FORMULA_H
