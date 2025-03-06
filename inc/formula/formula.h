#ifndef FORMULA_H

#include <cmath>  // 用于数学计算
#include <iostream>
#include <vector>

namespace formula
{

#include <algorithm>
#include <cmath>
#include <iostream>

const int N_S = 12;           // 滑移系统数量
const double k_B = 1.38e-23;  // 玻尔兹曼常数 (J/K)

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

int main()
{
    // 初始化物理参数
    double T = 300;          // 温度 (K)
    double P = 1e5;          // 压力 (Pa)
    double mu_0 = 26.1e9;    // 参考剪切模量 (Pa)
    double mu_p = 6.52e-2;   // 压力影响系数
    double mu_T = -6.16e-4;  // 温度影响系数
    double T_r = 300;        // 参考温度 (K)
    double mu = computeShearModulus(mu_0, mu_p, mu_T, P, T, T_r);

    double b = 2.86e-10;     // Burgers 矢量 (m)
    double v_g = 1.0e13;     // 位错尝试频率 (s^-1)
    double tau_ath = 1.0e7;  // 无热临界应力 (Pa)

    // 温度相关阻尼系数
    double B_T = 1.051e-5 * (1 + 0.9 * (T / 230));

    // 各滑移系统的剪应力 (Pa)
    double tau[N_S] = {-50e6, -45e6, -55e6, -40e6, -60e6, -52e6, -47e6, -53e6, -41e6, -59e6, -48e6, -51e6};

    // 滑移系统的初始位错密度 (m^-2)
    double rho_m[N_S] = {1.0e12,  0.9e12,  1.1e12, 1.2e12,  0.8e12, 1.05e12,
                         0.95e12, 1.15e12, 1.0e12, 1.05e12, 0.9e12, 1.1e12};
    double rho_im = 5.0e12;  // 非可动位错密度 (m^-2)
    double dt = 1e-6;        // 时间步长 (s)

    double alpha_het = 1.0e6;
    double alpha_mul = 0.1;
    double alpha_ann = 0.01;
    double alpha_trap = 0.02;
    double alpha_depin = 0.03;

    // 存储剪切应变率
    double gamma_dot[N_S];

    for (int i = 0; i < N_S; i++)
    {
        double delta_G = 0.5 * b * b * mu;  // 激活能量 (J)
        double t_w = computeDislocationWaitingTime(v_g, delta_G, T);
        double v_r = computeDislocationRunVelocity(tau[i], tau_ath, b, B_T);

        double L = 1.0 / sqrt(rho_m[i]);  // 计算位错间距 (m)
        double t_r = L / v_r;             // 计算位错运行时间 (s)
        double v = computeDislocationVelocity(L, t_w, t_r);

        // 计算剪切应变率
        gamma_dot[i] = computeShearStrainRate(rho_m[i], b, v, tau[i]);

        // 更新可动位错密度
        rho_m[i] = computeMobileDislocationDensity(rho_m[i], rho_im, gamma_dot[i], alpha_het, alpha_mul, alpha_ann,
                                                   alpha_trap, alpha_depin, dt);
    }

    // 输出结果
    std::cout << "剪切应变率 (s^-1) for each slip system:\n";
    for (int i = 0; i < N_S; i++)
    {
        std::cout << "γ_dot[" << i << "] = " << gamma_dot[i] << " s^-1\n";
    }

    return 0;
}

}  // namespace formula

#endif  // !FORMULA_H
