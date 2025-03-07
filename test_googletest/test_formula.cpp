#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace Eigen;

namespace slip_system
{

constexpr double k_B = 1.380649e-23;  // Boltzmann常数 (J/K)

// 材料参数（严格对应表1、表2）
struct MaterialParams
{
    // 热力学参数（表1）
    double C11 = 106.8e9, C12 = 60.4e9, C44 = 28.3e9;
    double rho0 = 2700.0;      // kg/m³
    double mu0 = 26.1e9;       // Pa
    double mu_p = 6.52e-2;     // 无量纲
    double mu_T = -6.16e-4;    // K⁻¹
    double cv = 900.0;         // J/(kg·K)
    double alpha_v = 69.9e-6;  // 热膨胀系数
    double beta_TQ = 0.9;      // Taylor-Quinney系数

    // 晶体塑性参数（表2）
    double tau_ath = 0.01e9;  // Pa
    double b = 2.86e-10;      // m
    double v_g = 1.0e13;      // s⁻¹
    double g0 = 0.5;          // 无量纲
    double p = 0.4, q = 2.5;
    double alpha_g = 0.7;
    double B = 1.051e-5;   // Pa·s
    double theta = 230.0;  // K
    double a_n[5] = {0.0, 1.0, 0.0, 0.0, 9.5e-3};
    double alpha_het = 1.2886e14;  // m⁻²
    double m_het = 0.85;
    double tau_a = 20.0e6;     // Pa
    double tau_b = 1.80668e9;  // Pa
    double alpha_mult = 0.116;
    double alpha_ann = 40.0;
    double alpha_trap = 21.0;    // m⁻¹
    double rho0_initial = 1e11;  // m⁻²
    double f_m = 0.5;            // 初始可动位错比例
    double m_rec = 1.0;
    double delta_s0 = 5e14;  // m⁻²
    double gamma0 = 1e7;     // s⁻¹
    double alpha_rec = 0.9;
};

// 晶体状态变量
struct CrystalState
{
    Matrix3d F_e = Matrix3d::Identity();
    Matrix3d F_p = Matrix3d::Identity();
    double T = 300.0;  // K
    double rho_m;      // 可动位错密度 (m⁻²)
    double rho_im;     // 不可动位错密度 (m⁻²)
    std::vector<double> gamma_dots;

    explicit CrystalState(const MaterialParams& mp)
        : rho_m(mp.rho0_initial * mp.f_m), rho_im(mp.rho0_initial * (1 - mp.f_m))
    {
    }
};

class CrystalPlasticityModel
{
   public:
    // 构造函数，使用MaterialParams参数初始化mp_和slip_systems_
    explicit CrystalPlasticityModel(const MaterialParams& mp) : mp_(mp), slip_systems_(generateFCCSlipSystems()) {}

    // 更新晶粒状态
    void updateState(CrystalState& state, const Matrix3d& F_total, double dt)
    {
        // 计算弹性应变增量
        Matrix3d F_e_trial = F_total * state.F_p.inverse();
        Matrix3d E_e = 0.5 * (F_e_trial.transpose() * F_e_trial - Matrix3d::Identity());
        // 计算第二PK应力
        Matrix3d S = computeSecondPKStress(E_e, state.T);
        // 计算解析切应力
        std::vector<double> tau_alpha = computeResolvedShearStresses(S, F_e_trial);
        // 计算滑移速率
        state.gamma_dots = computeSlipRates(tau_alpha, state);
        // 更新位错密度
        updateDislocationDensity(state, dt);
        // 更新温度
        updateTemperature(state, S, E_e, dt);
        // 更新塑性变形
        updatePlasticDeformation(state, dt);
    }

    // 计算第二PK应力
    Matrix3d computeSecondPKStress(const Matrix3d& E_e, double T)
    {
        // 定义应力矩阵
        Matrix3d S;
        // 计算E_e的迹
        double trE = E_e.trace();
        // 计算应力矩阵
        S = mp_.C11 * E_e + mp_.C12 * trE * Matrix3d::Identity();
        S += 0.5 * (mp_.C11 - mp_.C12) * (E_e * E_e).trace() * Matrix3d::Identity();
        // 计算温度对应力的影响
        S -= (mp_.alpha_v * (mp_.C11 + 2 * mp_.C12) / 3.0) * (T - 300.0) * Matrix3d::Identity();
        // 返回应力矩阵
        return S;
    }

   private:
    MaterialParams mp_;
    std::vector<std::pair<Vector3d, Vector3d>> slip_systems_;

    // 生成FCC滑移系统
    std::vector<std::pair<Vector3d, Vector3d>> generateFCCSlipSystems()
    {
        // 定义滑移系统
        std::vector<std::pair<Vector3d, Vector3d>> systems;
        // 定义法向量
        const std::vector<Vector3d> normals = {Vector3d(1, 1, 1).normalized(), Vector3d(1, -1, -1).normalized(),
                                               Vector3d(-1, 1, -1).normalized(), Vector3d(-1, -1, 1).normalized()};
        // 定义方向向量
        const std::vector<Vector3d> directions = {Vector3d(1, -1, 0).normalized(), Vector3d(1, 0, -1).normalized(),
                                                  Vector3d(0, 1, -1).normalized()};

        // 遍历法向量和方向向量，生成滑移系统
        for (const auto& n : normals)
        {
            for (const auto& s : directions)
            {
                systems.emplace_back(n, s);
                systems.emplace_back(n, Vector3d(s(1), s(2), s(0)));
                systems.emplace_back(n, Vector3d(s(2), s(0), s(1)));
            }
        }
        return systems;
    }

    // 计算解析剪切应力
    std::vector<double> computeResolvedShearStresses(const Matrix3d& S, const Matrix3d& F_e)
    {
        // 定义一个存储剪切应力的向量
        std::vector<double> tau_alpha;
        // 计算应力张量sigma
        Matrix3d sigma = (1.0 / F_e.determinant()) * F_e * S * F_e.transpose();

        // 遍历每个滑移系统
        for (const auto& [n, s] : slip_systems_)
        {
            // 计算Schmid因子
            Matrix3d schmid = 0.5 * (s * n.transpose() + n * s.transpose());
            // 计算剪切应力并添加到向量中
            tau_alpha.push_back(sigma.cwiseProduct(schmid).sum());
        }
        // 返回剪切应力向量
        return tau_alpha;
    }

    // 计算滑移速率
    std::vector<double> computeSlipRates(const std::vector<double>& tau_alpha, CrystalState& state)
    {
        // 定义滑移速率向量
        std::vector<double> gamma_dots;
        // 计算剪切模量
        const double mu = computeShearModulus(state);
        // 计算C_t
        const double C_t = sqrt(mu / mp_.rho0);

        // 遍历每个滑移系统
        for (size_t i = 0; i < slip_systems_.size(); ++i)
        {
            // 获取当前滑移系统的切应力
            const double tau = tau_alpha[i];
            // 计算g_alpha
            const double g_alpha = mp_.alpha_g * mu * mp_.b * sqrt(state.rho_m + state.rho_im);
            // 计算比例
            const double ratio = std::clamp((std::abs(tau) - mp_.tau_ath) / g_alpha, 0.0, 1.0);
            // 计算delta_G
            const double delta_G = mp_.g0 * mu * pow(mp_.b, 3) * pow(1 - pow(ratio, mp_.p), mp_.q);
            // 计算t_w
            const double t_w = 1.0 / mp_.v_g * (exp(delta_G / (k_B * state.T)) - 1.0);

            // 计算B_T
            double B_T = 0.0;
            for (int n = 0; n < 5; ++n)
                B_T += mp_.a_n[n] * pow(state.T / mp_.theta, n);
            B_T *= mp_.B;

            // 计算v_r和L
            double v_r = 0.0, L = 1.0 / sqrt(state.rho_im);
            // 迭代计算v_r
            for (int iter = 0; iter < 100; ++iter)
            {
                double v_new = (std::abs(tau) - mp_.tau_ath) * mp_.b / B_T * (1 - pow(v_r / C_t, 2));
                if (std::abs(v_new - v_r) < 1e-6)
                    break;
                v_r = v_new;
            }

            // 计算gamma_dot
            double gamma_dot = state.rho_m * mp_.b * L / (t_w + L / (v_r + 1e-20));
            gamma_dot *= (tau >= 0 ? 1.0 : -1.0);
            gamma_dots.push_back(gamma_dot);
        }
        return gamma_dots;
    }

    // 更新晶体的位错密度
    void updateDislocationDensity(CrystalState& state, double dt)
    {
        // 计算总位错密度
        double total_gamma = 0.0;
        for (double g : state.gamma_dots)
            total_gamma += std::abs(g) * dt;

        // 计算异质位错密度
        double rho_het = 0.0;
        if (total_gamma >= mp_.tau_a && total_gamma <= mp_.tau_b)
        {
            double x = (total_gamma - mp_.tau_a) / (mp_.tau_b - mp_.tau_a);
            rho_het = mp_.alpha_het * (mp_.m_het + 1) * pow(x, mp_.m_het);
        }

        // 计算多重位错密度
        double rho_mult = mp_.alpha_mult * total_gamma / (computeShearModulus(state) * pow(mp_.b, 2));
        // 计算湮灭位错密度
        double rho_ann = mp_.alpha_ann * (state.rho_m + state.rho_im) * total_gamma;
        // 计算俘获位错密度
        double rho_trap = mp_.alpha_trap * sqrt(state.rho_m) * total_gamma;

        // 计算饱和度
        double delta_sat = mp_.delta_s0 * pow(total_gamma / mp_.gamma0, mp_.alpha_rec);
        // 计算脱捕获位错密度
        double alpha_depin =
            pow(state.rho_im / sqrt(pow(delta_sat, 2) + pow(mp_.delta_s0, 2)), mp_.m_rec) * mp_.alpha_trap;
        double rho_depin = alpha_depin * sqrt(state.rho_im) * total_gamma;

        // 更新位错密度
        state.rho_m += rho_het + rho_mult - rho_ann - rho_trap + rho_depin;
        state.rho_im += rho_trap - rho_depin;
        // 确保位错密度非负
        state.rho_m = std::max(state.rho_m, 0.0);
        state.rho_im = std::max(state.rho_im, 0.0);
    }

    // 更新晶体的温度
    void updateTemperature(CrystalState& state, const Matrix3d& S, const Matrix3d& E_e, double dt)
    {
        // 初始化塑性应变能
        double W_p = 0.0;
        // 遍历所有滑移系统
        for (size_t i = 0; i < slip_systems_.size(); ++i)
        {
            // 获取滑移系统的法向量和滑移方向
            const auto& [n, s] = slip_systems_[i];
            // 计算塑性应变能
            W_p += S.cwiseProduct(s * n.transpose()).sum() * state.gamma_dots[i] * dt;
        }
        // 计算温度变化的第一项
        double term1 = mp_.beta_TQ * W_p / (mp_.rho0 * mp_.cv);
        // 计算温度变化的第二项
        double term2 = state.T * mp_.alpha_v * E_e.trace() / (mp_.rho0 * mp_.cv);
        // 更新温度
        state.T += term1 - term2;
    }

    // 更新晶体的塑性变形
    void updatePlasticDeformation(CrystalState& state, double dt)
    {
        // 初始化塑性变形矩阵L_p为0矩阵
        Matrix3d L_p = Matrix3d::Zero();
        // 遍历所有滑移系统
        for (size_t i = 0; i < slip_systems_.size(); ++i)
        {
            // 获取滑移系统的法向量和滑移矢量
            const auto& [n, s] = slip_systems_[i];
            // 计算塑性变形矩阵L_p
            L_p += state.gamma_dots[i] * s * n.transpose();
        }
        // 更新塑性变形矩阵F_p
        state.F_p = (Matrix3d::Identity() + dt * L_p) * state.F_p;
    }

    // 计算切变模量
    double computeShearModulus(const CrystalState& state)
    {
        // 计算应变矩阵的行列式
        double J_e = state.F_e.determinant();
        // 返回切变模量，其中mp_为材料参数，mu0为切变模量，mu_p为应变引起的切变模量，mu_T为温度引起的切变模量
        return mp_.mu0 * (1.0 + mp_.mu_p * J_e + mp_.mu_T * (state.T - 300.0));
    }
};

}  // namespace slip_system

TEST(CrystalPlasticity, StressOvershoot)
{
    using namespace slip_system;

    MaterialParams mp;
    CrystalState state(mp);
    CrystalPlasticityModel model(mp);

    Matrix3d F = Matrix3d::Identity();
    const double strain_rate = 1e7;
    const double total_strain = 0.1;
    const int steps = 100;
    const double dt = total_strain / (strain_rate * steps);

    std::vector<double> stress_history;
    for (int i = 0; i < steps; ++i)
    {
        F(2, 2) = 1.0 - strain_rate * dt * (i + 1);
        model.updateState(state, F, dt);

        Matrix3d E_e = 0.5 * (state.F_e.transpose() * state.F_e - Matrix3d::Identity());
        Matrix3d S = model.computeSecondPKStress(E_e, state.T);

        // 正确计算冯·米塞斯应力
        double s11_s22 = S(0, 0) - S(1, 1);
        double s22_s33 = S(1, 1) - S(2, 2);
        double s33_s11 = S(2, 2) - S(0, 0);
        double von_mises = sqrt(0.5 * (s11_s22 * s11_s22 + s22_s33 * s22_s33 + s33_s11 * s33_s11));

        stress_history.push_back(von_mises);
    }

    // 验证应力超调现象
    auto max_it = std::max_element(stress_history.begin(), stress_history.end());
    EXPECT_GT(*max_it, stress_history.back())
        << "未检测到应力超调，最大应力：" << *max_it << "，最终应力：" << stress_history.back();
}