以下是论文中所有公式与代码实现的详细对应关系，按出现顺序整理：

---

### **1. 运动学分解**
**公式 (1.1)**
$\mathbf{F} = \mathbf{F}^e \mathbf{F}^p$
```cpp
// 在updateState()中体现
Matrix3d F_e_trial = F_total * state.F_p.inverse(); 
```

---

### **2. 塑性速度梯度**
**公式 (1.2)**
$\mathbf{L}^p = \sum_{\alpha} \dot{\gamma}^{\alpha} \mathbf{s}_0^{\alpha} \otimes \mathbf{n}_0^{\alpha}$
```cpp
// updatePlasticDeformation()中
for (size_t alpha...)
    L_p += gamma_dots[alpha] * (s * n.transpose());
```

---

### **3. 热力学第二定律**
**公式 (1.3)-(1.8)**
Clausius-Duhem不等式的推导
```cpp
// 隐含在应力计算和能量守恒中，无直接对应代码
```

---

### **4. Helmholtz自由能展开**
**公式 (1.9)**
$\rho \varphi = \frac{1}{2}C_{ijkl}E^e_{ij}E^e_{kl} + ...$
```cpp
Matrix3d computeSecondPKStress(...) {
    // 二阶项
    S = C11*E_e + C12*trace... 
    // 三阶项
    S(i,j) += 0.5*C111*E_e(i,j)^3 + ... 
}
```

---

### **5. 二阶PK应力**
**公式 (1.10)**
$\mathbf{S} = \rho \frac{\partial \varphi}{\partial \mathbf{E}^e}$
```cpp
Matrix3d computeSecondPKStress(...) {
    // 完整实现三阶弹性项
}
```

---

### **6. 温度演化方程**
**公式 (1.13)**
$\dot{T} = \frac{\beta_T}{\rho c_V} \mathbf{S}:\dot{\mathbf{E}}^p - \frac{T \alpha:\dot{\mathbf{E}}^e}{\rho c_V}$
```cpp
void updateTemperature(...) {
    term1 = beta_TQ * W_p / (rho0*cv);
    term2 = T*alpha_v*(E_e.trace())/...
    T += (term1 - term2)*dt;
}
```

---

### **7. Orowan方程**
**公式 (1.14)**
$\dot{\gamma}^{\alpha} = \rho_m^{\alpha} b \bar{v}^{\alpha} \text{sgn}(\tau^{\alpha})$
```cpp
// computeSlipRates()中
gamma_dot = rho_m * b * L / (t_w + L/v_r);
```

---

### **8. 位错速度模型**
**公式 (1.15)-(1.20)**
$\bar{v}^{\alpha} = \frac{L^{\alpha}}{t_w^{\alpha} + t_r^{\alpha}}$
```cpp
// computeSlipRates()中的迭代计算
for (iter...) {
    v_new = (|tau|-tau_ath)*b / (B_T*(1+(v_r/C_t)^2);
}
```

---

### **9. 位错密度演化**
**公式 (1.21)-(1.30)**
$\dot{\rho}_m^{\alpha} = \dot{\rho}_{nuc}^{\alpha} + ...$
```cpp
void updateDislocationDensity(...) {
    // 异质成核
    if (tau in [tau_a,tau_b])...
    // 增殖、湮灭等项
    rho_mult = alpha_mult * total_gamma / b²...
}
```

---

### **10. 塑性变形梯度更新**
**公式 (1.32)-(1.34)**
$\tilde{\mathbf{F}}_{t+\Delta t}^{p-1} = \tilde{\mathbf{F}}_t^{p-1}(I - \sum \Delta \gamma^{\alpha} \mathbf{S}_0^{\alpha})$
```cpp
void updatePlasticDeformation(...) {
    // SVD分解保持det(F_p)=1
    JacobiSVD<Matrix3d> svd(...);
    S = normalized(S); 
}
```

---

### **未直接实现的公式**
1. **公式 (1.4)-(1.8)**：热力学不等式推导，仅通过能量守恒间接体现
2. **公式 (1.11)-(1.12)**：热传导方程，因假设绝热过程未显式实现
3. **公式 (1.31)**：晶格坐标系变换，隐含在滑移系定义中

---

### **关键参数对照**
| 论文参数        | 代码变量         | 位置               |
|----------------|----------------|-------------------|
| $\mu_0$        | mp_.mu0        | MaterialParams    |
| $B^*(T)$       | B_T            | computeSlipRates()|
| $\alpha_{mult}$| mp_.alpha_mult | updateDislocationDensity()|

---

### **需要特别注意的实现差异**
1. **三阶弹性项简化**：代码中的三阶项实现采用了张量分量循环展开，而非论文中的完整爱因斯坦求和约定
2. **压力计算近似**：`computeShearModulus()`中使用$\text{tr}(S)/3$近似压力，而非严格通过状态方程
3. **位错速度迭代**：代码设置最大迭代次数50次，而论文未明确说明迭代方法

建议在关键物理量（如位错密度、温度）添加更多验证测试，确保与论文结果的一致性。