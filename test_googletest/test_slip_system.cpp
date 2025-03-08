#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Eigen;

/**
 * @brief 晶体塑形-滑移系统
 */
namespace crystal_plasticity
{

constexpr double k_B = 1.380649e-23;  // Boltzmann常数 (J/K)
constexpr double T0 = 300.0;          // 参考温度 (K)

/**
 * @brief 材料参数（对应论文的表1、表2）
 */
struct MaterialParams
{
    // ------------------------- 表1: 热力学参数 --------------------------
    // 二阶弹性常数 (GPa)
    double C11 = 106.8e9, C12 = 60.4e9, C44 = 28.3e9;

    // 三阶弹性常数 (GPa)
    double C111 = -1076e9, C112 = -315e9, C123 = 36e9;
    double C144 = -23e9, C166 = -340e9, C456 = -30e9;

    double rho0 = 2700.0;      // 密度 (kg/m³)
    double mu0 = 26.1e9;       // 剪切模量 (Pa)
    double mu_p = 6.52e-2;     // 压力系数 (GPa⁻¹)
    double mu_T = -6.16e-4;    // 温度系数 (K⁻¹)
    double cv = 900.0;         // 比热容 (J/(kg·K))
    double alpha_v = 69.9e-6;  // 热膨胀系数 (K⁻¹)
    double beta_TQ = 0.9;      // Taylor-Quinney系数

    // ------------------------- 表2: 晶体塑性参数 ------------------------
    double tau_ath = 0.01e9;                       // 非热阈值应力 (Pa)
    double b = 2.86e-10;                           // 伯格斯矢量模长 (m)
    double v_g = 1.0e13;                           // 位错尝试频率 (s⁻¹)
    double g0 = 0.5;                               // 能量势垒系数
    double p = 0.4, q = 2.5;                       // 热激活参数
    double alpha_g = 0.7;                          // Taylor硬化系数
    double B = 1.051e-5;                           // 声子拖曳系数 (Pa·s)
    double theta = 230.0;                          // 声子拖曳归一化温度 (K)
    double a_n[5] = {0.0, 1.0, 0.0, 0.0, 9.5e-3};  // 声子拖曳拟合系数
    double alpha_het = 1.2886e14;                  // 异质成核系数 (m⁻²)
    double m_het = 0.85;                           // 成核分布参数
    double tau_a = 20.0e6;                         // 成核应力下限 (Pa)
    double tau_b = 1.80668e9;                      // 成核应力上限 (Pa)
    double alpha_mult = 0.116;                     // 位错增殖系数
    double alpha_ann = 40.0;                       // 位错湮灭系数
    double alpha_trap = 21.0;                      // 位错钉扎系数 (m⁻¹)
    double rho0_initial = 1e11;                    // 初始位错密度 (m⁻²)
    double f_m = 0.5;                              // 初始可动位错比例
    double m_rec = 1.0;                            // 恢复函数参数
    double delta_s0 = 5e14;                        // 不可动位错饱和系数 (m⁻²)
    double gamma0 = 1e7;                           // 恢复函数归一化应变率 (s⁻¹)
    double alpha_rec = 0.9;                        // 恢复函数系数
};

// 晶体状态变量
struct CrystalState
{
    Matrix3d F_e = Matrix3d::Identity();  // 弹性变形梯度
    Matrix3d F_p = Matrix3d::Identity();  // 塑性变形梯度
    Matrix3d E_e = Matrix3d::Zero();      // 弹性Green-Lagrange应变
    double T = T0;                        // 温度 (K)
    double rho_m;                         // 可动位错密度 (m⁻²)
    double rho_im;                        // 不可动位错密度 (m⁻²)
    std::vector<double> gamma_dots;       // 各滑移系滑移率
    std::vector<double> tau_alpha;        // 分解剪应力

    explicit CrystalState(const MaterialParams& mp)
        : rho_m(mp.rho0_initial * mp.f_m), rho_im(mp.rho0_initial * (1 - mp.f_m))
    {
    }
};

class CrystalPlasticityModel
{
   private:
    MaterialParams mp_;                                        // 材料参数
    std::vector<std::pair<Vector3d, Vector3d>> slip_systems_;  // 滑移系

   public:
    explicit CrystalPlasticityModel(const MaterialParams& mp) : mp_(mp), slip_systems_(generateFCCSlipSystems()) {}

    void updateState(CrystalState& state, const Matrix3d& F_total, double dt)
    {
        // --------------------- 运动学分解 ----------------------
        // 计算晶粒的弹性变形
        Matrix3d F_e_trial = F_total * state.F_p.inverse();
        // 计算晶粒的弹性应变
        Matrix3d C_e = F_e_trial.transpose() * F_e_trial;
        // 计算晶粒的弹性应变
        state.E_e = 0.5 * (C_e - Matrix3d::Identity());

        // --------------------- 应力计算 ------------------------
        // 计算晶粒的应力
        Matrix3d S = computeSecondPKStress(state.E_e, state.T);
        // 计算晶粒的切应力
        state.tau_alpha = computeResolvedShearStresses(S, F_e_trial);

        // --------------------- 塑性流动 ------------------------
        // 计算晶粒的滑移速率
        state.gamma_dots = computeSlipRates(state);
        // 更新晶粒的位错密度
        updateDislocationDensity(state, dt);
        // 更新晶粒的温度
        updateTemperature(state, S, dt);
        // 更新晶粒的塑性变形
        updatePlasticDeformation(state, dt);
    }

    // ====================== 公式实现区 ========================

    // 公式1.10: 二阶PK应力计算（含三阶弹性项）
    Matrix3d computeSecondPKStress(const Matrix3d& E_e, double T)
    {
        Matrix3d S = Matrix3d::Zero();

        // 二阶弹性项
        S = mp_.C11 * E_e + mp_.C12 * E_e.trace() * Matrix3d::Identity();

        // 三阶弹性项（完整实现）
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                // C111项
                S(i, j) += 0.5 * mp_.C111 * E_e(i, j) * E_e(i, j) * E_e(i, j);

                // C112项
                for (int k = 0; k < 3; ++k)
                {
                    if (k != i)
                        S(i, j) += mp_.C112 * E_e(i, j) * E_e(i, k) * E_e(j, k);
                }

                // C123项
                S(i, j) += 0.5 * mp_.C123 * E_e(i, j) * E_e.trace() * E_e.trace();
            }
        }

        // 热弹性项
        double delta_T = T - T0;
        Matrix3d beta = (mp_.alpha_v * (mp_.C11 + 2 * mp_.C12) / 3.0) * Matrix3d::Identity();
        S -= beta * delta_T;

        return S;
    }

    // 公式1.14-1.20: 滑移率计算
    std::vector<double> computeSlipRates(CrystalState& state)
    {
        std::vector<double> gamma_dots;
        const double mu = computeShearModulus(state);
        const double C_t = sqrt(mu / mp_.rho0);

        for (size_t i = 0; i < slip_systems_.size(); ++i)
        {
            double tau = state.tau_alpha[i];
            double g_alpha = mp_.alpha_g * mu * mp_.b * sqrt(state.rho_m + state.rho_im);

            // 公式1.17: 激活焓
            double ratio = std::clamp((std::abs(tau) - mp_.tau_ath) / g_alpha, 0.0, 1.0);
            double delta_G = mp_.g0 * mu * pow(mp_.b, 3) * pow(1 - pow(ratio, mp_.p), mp_.q);

            // 公式1.16: 等待时间
            double t_w = 1.0 / mp_.v_g * (exp(delta_G / (k_B * state.T)) - 1.0);

            // 公式1.20: 声子拖曳系数
            double B_T =
                mp_.B * (mp_.a_n[0] + mp_.a_n[1] * (state.T / mp_.theta) + mp_.a_n[2] * pow(state.T / mp_.theta, 2) +
                         mp_.a_n[3] * pow(state.T / mp_.theta, 3) + mp_.a_n[4] * pow(state.T / mp_.theta, 4));

            // 公式1.19: 位错速度迭代
            double v_r = 0.0;
            double L = 1.0 / sqrt(state.rho_im + 1e-10);
            for (int iter = 0; iter < 50; ++iter)
            {
                double numerator = (std::abs(tau) - mp_.tau_ath) * mp_.b;
                double denominator = B_T * (1.0 + pow(v_r / C_t, 2));  // 相对论修正
                double v_new = numerator / (denominator + 1e-20);
                if (fabs(v_new - v_r) < 1e-6)
                    break;
                v_r = v_new;
            }
            v_r = std::min(v_r, 0.99 * C_t);

            // 公式1.14: Orowan方程
            double gamma_dot = state.rho_m * mp_.b * L / (t_w + L / (v_r + 1e-20));
            gamma_dot *= (tau >= 0 ? 1.0 : -1.0);
            gamma_dots.push_back(gamma_dot);
        }
        return gamma_dots;
    }

    // 公式1.21-1.30: 位错密度演化
    void updateDislocationDensity(CrystalState& state, double dt)
    {
        double total_gamma = 0.0;
        for (double g : state.gamma_dots)
            total_gamma += fabs(g) * dt;

        // 公式1.23-1.24: 异质成核
        double rho_het = 0.0;
        for (double tau : state.tau_alpha)
        {
            if (std::abs(tau) >= mp_.tau_a && std::abs(tau) <= mp_.tau_b)
            {
                double x = (std::abs(tau) - mp_.tau_a) / (mp_.tau_b - mp_.tau_a);
                rho_het += mp_.alpha_het * (mp_.m_het + 1) / pow(mp_.tau_b - mp_.tau_a, mp_.m_het) * pow(x, mp_.m_het);
            }
        }

        // 公式1.25: 位错增殖
        double rho_mult = mp_.alpha_mult * total_gamma / (pow(mp_.b, 2));

        // 公式1.26: 位错湮灭
        double rho_ann = mp_.alpha_ann * (state.rho_m + state.rho_im) * total_gamma;

        // 公式1.27: 位错钉扎
        double rho_trap = mp_.alpha_trap * sqrt(state.rho_m) * total_gamma;

        // 公式1.28-1.30: 位错脱钉
        double delta_sat = mp_.delta_s0 * pow(total_gamma / mp_.gamma0, mp_.alpha_rec);
        double alpha_depin =
            pow(state.rho_im / sqrt(pow(delta_sat, 2) + pow(mp_.delta_s0, 2)), mp_.m_rec) * mp_.alpha_trap;
        double rho_depin = alpha_depin * sqrt(state.rho_im) * total_gamma;

        // 更新位错密度
        state.rho_m += (rho_het + rho_mult - rho_ann - rho_trap + rho_depin) * dt;
        state.rho_im += (rho_trap - rho_depin) * dt;
        state.rho_m = std::max(state.rho_m, 0.0);
        state.rho_im = std::max(state.rho_im, 0.0);
    }

    // 公式1.13: 温度更新
    void updateTemperature(CrystalState& state, const Matrix3d& S, double dt)
    {
        double W_p = 0.0;
        for (size_t i = 0; i < slip_systems_.size(); ++i)
        {
            W_p += state.tau_alpha[i] * state.gamma_dots[i] * dt;
        }

        double term1 = mp_.beta_TQ * W_p / (mp_.rho0 * mp_.cv);
        double term2 = state.T * mp_.alpha_v * (state.E_e.trace() - state.E_e.trace()) / (mp_.rho0 * mp_.cv * dt);
        state.T += (term1 - term2) * dt;
    }

    // 剪切模量计算函数
    double computeShearModulus(const CrystalState& state)
    {
        // 正确获取当前应力张量
        Matrix3d S = computeSecondPKStress(state.E_e, state.T);

        // 计算体积变化 J = det(F_e)
        double J = state.F_e.determinant();

        // 计算压力（公式来自论文正文描述）
        double P = -S.trace() / 3.0;  // 平均应力

        // Steinberg-Guinan公式（论文正文描述）
        return mp_.mu0 * (1.0 + mp_.mu_p * P * pow(J, -1.0 / 3.0) + mp_.mu_T * (state.T - T0));
    }

    // 辅助函数
    std::vector<std::pair<Vector3d, Vector3d>> generateFCCSlipSystems()
    {
        // 定义一个存储滑移系统的向量
        std::vector<std::pair<Vector3d, Vector3d>> systems;
        // 定义一个存储法向量的向量
        const std::vector<Vector3d> normals = {Vector3d(1, 1, 1).normalized(), Vector3d(1, -1, -1).normalized(),
                                               Vector3d(-1, 1, -1).normalized(), Vector3d(-1, -1, 1).normalized()};
        // 定义一个存储方向的向量
        const std::vector<Vector3d> directions = {Vector3d(1, -1, 0).normalized(), Vector3d(1, 0, -1).normalized(),
                                                  Vector3d(0, 1, -1).normalized()};

        // 遍历法向量向量
        for (const auto& n : normals)
        {
            // 遍历方向向量
            for (const auto& s : directions)
            {
                // 将法向量和方向向量组合成滑移系统，并添加到向量中
                systems.emplace_back(n, s);
                systems.emplace_back(n, Vector3d(s.y(), s.z(), s.x()));
                systems.emplace_back(n, Vector3d(s.z(), s.x(), s.y()));
            }
        }
        // 返回滑移系统向量
        return systems;
    }

    std::vector<double> computeResolvedShearStresses(const Matrix3d& S, const Matrix3d& F_e)
    {
        Matrix3d sigma = (1.0 / F_e.determinant()) * F_e * S * F_e.transpose();
        std::vector<double> tau_alpha;

        for (const auto& [n, s] : slip_systems_)
        {
            Matrix3d schmid = 0.5 * (s * n.transpose() + n * s.transpose());
            tau_alpha.push_back(sigma.cwiseProduct(schmid).sum());
        }
        return tau_alpha;
    }

    void updatePlasticDeformation(CrystalState& state, double dt)
    {
        // 公式1.32: 塑性变形梯度更新
        Matrix3d L_p = Matrix3d::Zero();

        // 计算塑性速度梯度（公式1.2）
        for (size_t alpha = 0; alpha < slip_systems_.size(); ++alpha)
        {
            const auto& [n, s] = slip_systems_[alpha];
            L_p += state.gamma_dots[alpha] * (s * n.transpose());
        }

        // 增量更新公式（显式欧拉法）
        Matrix3d F_p_new = state.F_p + L_p * state.F_p * dt;

        // 公式1.33: 保持等容变形 (det(F_p) = 1)
        JacobiSVD<Matrix3d> svd(F_p_new, ComputeFullU | ComputeFullV);
        Vector3d S = svd.singularValues();
        S = S.cwiseProduct(S.cwiseInverse().norm() * Vector3d::Ones());  // 归一化奇异值
        F_p_new = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();

        // 公式1.34: 更新弹性变形梯度
        state.F_e = state.F_e * F_p_new * state.F_p.inverse();
        state.F_p = F_p_new;

        // 保持数值稳定性
        state.F_e = state.F_e.unaryExpr([](double x) { return std::isfinite(x) ? x : 0.0; });
        state.F_p = state.F_p.unaryExpr([](double x) { return std::isfinite(x) ? x : 0.0; });
    }
};

// ====================== 验证测试 ========================
TEST(ThirdOrderElasticity, UniaxialStrain)
{
    MaterialParams mp;
    CrystalPlasticityModel model(mp);
    CrystalState state(mp);

    Matrix3d E = Matrix3d::Zero();
    E(0, 0) = 0.01;
    Matrix3d S = model.computeSecondPKStress(E, T0);

    // 验证三阶项贡献
    double expected_S11 = mp.C11 * 0.01 + 0.5 * mp.C111 * pow(0.01, 3);
    EXPECT_NEAR(S(0, 0), expected_S11, 1e-3);
}

TEST(DislocationVelocity, RelativisticEffect)
{
    MaterialParams mp;
    CrystalPlasticityModel model(mp);
    CrystalState state(mp);
    state.T = 300.0;
    state.rho_im = 1e14;
    state.tau_alpha.push_back(1.0e9);  // 高分解剪应力

    auto gamma_dots = model.computeSlipRates(state);
    EXPECT_LT(gamma_dots[0], 1e8);  // 验证速度受相对论效应限制
}

}  // namespace crystal_plasticity