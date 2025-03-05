#include <gtest/gtest.h>

#include "formula/formula.h"

using namespace formula;

TEST(FormulaTest, test)
{
    // 初始化物理参数
    double T = 300;          // 温度 (K)
    double P = 1e5;          // 压力 (Pa)
    double mu_0 = 26.1e9;    // 参考剪切模量 (Pa)
    double mu_p = 6.52e-2;   // 压力影响系数
    double mu_T = -6.16e-4;  // 温度影响系数
    double T_r = 300;        // 参考温度 (K)
    double mu = computeShearModulus(mu_0, mu_p, mu_T, P, T, T_r);

    double b = 2.86e-10;                // Burgers 矢量 (m)
    double tau = -50e6;                 // 剪应力 (Pa)
    double v_g = 1.0e13;                // 位错尝试频率 (s^-1)
    double delta_G = 0.5 * b * b * mu;  // 激活能量 (J)
    double t_w = computeDislocationWaitingTime(v_g, delta_G, T);

    double B_T = 1.051e-5 * (1 + 0.9 * (T / 230));  // 温度相关阻尼系数
    double tau_ath = 1.0e7;                         // 无热临界应力 (Pa)
    double v_r = computeDislocationRunVelocity(tau, tau_ath, b, B_T);

    // 计算可动位错密度 (初始值)
    double rho_m = 1.0e12;   // 可动位错密度 (m^-2)
    double rho_im = 5.0e12;  // 非可动位错密度 (m^-2)
    double dt = 1e-6;        // 时间步长 (s)
    double alpha_het = 1.0e6;
    double alpha_mul = 0.1;
    double alpha_ann = 0.01;
    double alpha_trap = 0.02;
    double alpha_depin = 0.03;

    double L = 1.0 / sqrt(rho_m);  // 计算位错间距 (m)
    double t_r = L / v_r;          // 计算位错运行时间 (s)
    double v = computeDislocationVelocity(L, t_w, t_r);

    // 计算剪切应变率 (s^-1)
    double gamma_dot = computeShearStrainRate(rho_m, b, v, tau);

    // 更新可动位错密度
    rho_m = computeMobileDislocationDensity(rho_m, rho_im, gamma_dot, alpha_het, alpha_mul, alpha_ann, alpha_trap,
                                            alpha_depin, dt);

    // 输出计算结果
    std::cout << "剪切应变率: " << gamma_dot << " s^-1" << std::endl;
    std::cout << "更新后的可动位错密度: " << rho_m << " m^-2" << std::endl;
}