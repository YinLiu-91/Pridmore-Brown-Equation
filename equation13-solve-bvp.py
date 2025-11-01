import numpy as np
from scipy.integrate import solve_bvp
from scipy.special import jv, jvp
import matplotlib.pyplot as plt


# ------------------------------------------------------------------------------------
# 参数设置（可自行修改）
# ------------------------------------------------------------------------------------
m = 24                  # 模态阶数
k = 31.0                # 无量纲波数
M = 0.5                 # 马赫数
kmn = 25.42 + 0.13j     # 轴向波数（可取复数）
Z = 3.0 - 0.5j          # 壁面阻抗（复数）
r_min = 1.0e-4          # 避免 r=0 的奇点，从一个很小的半径开始
r_max = 1.0


# ------------------------------------------------------------------------------------
# 将复数问题拆为实数 4 维系统
# y = [Re(p), Im(p), Re(dp/dr), Im(dp/dr)]
# ------------------------------------------------------------------------------------

def ode(r, y):
    p = y[0] + 1j * y[1]
    dp = y[2] + 1j * y[3]

    coeff = (k - M * kmn) ** 2 - kmn ** 2
    d2p = -(1.0 / r) * dp - (coeff - m ** 2 / r ** 2) * p

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
plt.xlabel('r')
plt.ylabel('声压')
plt.title(f'公式(13) BVP 解，m={m}')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()