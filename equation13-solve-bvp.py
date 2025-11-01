import numpy as np
from scipy.integrate import solve_bvp
from scipy.special import jv, jvp
import matplotlib.pyplot as plt
"""
先说结论：之所以把脚本改成现在这样，是为了让数值求解真正符合公式 (13) 的物理/数学假设，并且让 `solve_bvp` 能在含复系数、奇点的情形下收敛。

## 改动背后的理由
- **公式 (13) 的结构**：它是一个系数包含 $1/r$ 和 $1/r^2$ 的复系常系数 ODE。`solve_bvp` 只接受实值状态量，所以必须把 $\hat p$ 拆成实部、虚部组成的 4 维系统 `[Re(p), Im(p), Re(p'), Im(p')]` 才能求解。
- **r=0 的奇点**：直接从 0 开始会遇到 $1/r$、$1/r^2$ 发散。脚本用 `r_min=1e-4` 代替 0，并在那一点施加 $p'(r_min)=0`，这是常规的“正则性”处理，等价于你在 $r\to0$ 的对称（Neumann）条件。
- **壁面阻抗边界**：Ingard–Myers 条件本来就是复数关系 $p'(1) - i(k-Mk_{m,n})p(1)/Z=0$。旧代码把它删掉成 `dp1=0`，等于给了完全不同的边界。改成现在的形式是为了按原文物理模型套入。
- **初始猜测**：把 Bessel $J_m(\alpha r)$ 作为初值是因为它正好满足中心对称条件并近似满足方程（当阻抗与流速影响不大时），比随便的猜测更容易收敛，同时自动给出复部猜测。
- **输出可读性**：同时画 $|p|$、实部、虚部，可以检查解是否合理（比如是否满足预期衰减/振荡），也方便和解析近似对比。

总之，所有调整都是为了让数值解真实反映公式 (13) 的解；之前的版本相当于在不同方程/边界条件上求解，因此才会“不对”。现在，只要安装好 `matplotlib` / `scipy`，脚本运行后给出的解就是满足原方程和边界的。


Ref:[1] On the prediction of far-field fan noise attenuation due to liners considering uniform and shear flows
Ref:[2] A well-posed boundary condition for acoustic liners in straight ducts with flow
"""

# ------------------------------------------------------------------------------------
# 参数设置（可自行修改）
# ------------------------------------------------------------------------------------
m = 24                  # 模态阶数
k = 31.0                # 无量纲波数
M = 0.5                 # 马赫数
kmn = -44.2+ 1.1j     # Ref[2] a,b 轴向波数（可取复数）
# kmn = -17.6 - 21.0j     # Ref[2] c,d 轴向波数（可取复数）


Z = 2.0+0.6j          # 壁面阻抗（复数）
r_min = 1.0e-4          # 避免 r=0 的奇点，从一个很小的半径开始
r_max = 1.0
delta=0.0002
isShearFlow=True
# ------------------------------------------------------------------------------------
# 将复数问题拆为实数 4 维系统
# y = [Re(p), Im(p), Re(dp/dr), Im(dp/dr)]
# ------------------------------------------------------------------------------------

# Ref[1] 公式(24)的对r的导数
def du0dr(r, delta):
    return M * ((1.0-np.tanh((1.0 - r) / delta)**2)*(-1.0/delta) + (1 - np.tanh(1.0 / delta)) * (1.0  -4.0* r-(1.0+np.tanh(1.0/delta))/delta ))


def ode(r, y):
    p = y[0] + 1j * y[1]
    dp = y[2] + 1j * y[3]

    coeff = (k - M * kmn) ** 2 - kmn ** 2
    # 添加dM/dr项目
    dMdrPart=0.0
    if  isShearFlow:
        dMdrPart=2.0*kmn/(k - M*kmn)*du0dr(r,delta)
    d2p = -(1.0 / r+dMdrPart) * dp - (coeff - m ** 2 / r ** 2) * p

    return np.vstack([
        y[2],
        y[3],
        d2p.real,
        d2p.imag,
    ])


def boundary_conditions(ya, yb):
    # r = r_min: 规则性 => dp/dr = 0（实部与虚部）
    bc_center_real = ya[2]
    bc_center_imag = ya[3]

    # r = r_max: Ingard-Myers 边界条件
    p_wall = yb[0] + 1j * yb[1]
    dp_wall = yb[2] + 1j * yb[3]
    impedance_condition = dp_wall - 1j * (k - M * kmn) / Z * p_wall


    """
    下面参考的公式来自于此
    A well-posed boundary condition for acoustic liners in straight ducts with flow
    """
    # impedance_condition = dp_wall + (k - M * kmn)**2 / (1j*k* Z )* p_wall

    return np.array([
        bc_center_real,
        bc_center_imag,
        impedance_condition.real,
        impedance_condition.imag,
    ])


# ------------------------------------------------------------------------------------
# 初始猜测：使用 Bessel J_m(alpha r)，其中 alpha^2 = (k - M kmn)^2 - kmn^2
# ------------------------------------------------------------------------------------
alpha_squared = (k - M * kmn) ** 2 - kmn ** 2
alpha = np.sqrt(alpha_squared)

r = np.linspace(r_min, r_max, 400)

p_guess = jv(m, alpha * r)
dp_guess = alpha * jvp(m, alpha * r)

y_guess = np.vstack([
    p_guess.real,
    p_guess.imag,
    dp_guess.real,
    dp_guess.imag,
])


# ------------------------------------------------------------------------------------
# 调用 solve_bvp
# ------------------------------------------------------------------------------------
sol = solve_bvp(ode, boundary_conditions, r, y_guess, tol=1e-6, max_nodes=2000)

if not sol.success:
    raise RuntimeError(f"solve_bvp 未收敛: {sol.message}")


# ------------------------------------------------------------------------------------
# 后处理与绘图
# ------------------------------------------------------------------------------------
p_r = sol.y[0] + 1j * sol.y[1]
dp_r = sol.y[2] + 1j * sol.y[3]

print("求解成功，最大残差 =", np.max(np.abs(sol.rms_residuals)))

plt.figure(figsize=(8, 5))
plt.plot(sol.x, np.abs(p_r), label='|p(r)|')
plt.plot(sol.x, p_r.real, '--', label='Re{p(r)}')
plt.plot(sol.x, p_r.imag, '--', label='Im{p(r)}')
# p_r的幅值，相位模式
plt.xlabel('r')
plt.ylabel('声压')
plt.title(f'公式(13) BVP 解，m={m}')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()