import numpy as np
from scipy.integrate import solve_bvp
from scipy.special import jv, jvp
import matplotlib.pyplot as plt

"""
Ref:[1] On the prediction of far-field fan noise attenuation due to liners considering uniform and shear flows
Ref:[2] A well-posed boundary condition for acoustic liners in straight ducts with flow
1. 主要参考了Ref[1]
2. 根据Ref[2]做了修改
3. 结果与Ref[2]做对比
"""

# ------------------------------------------------------------------------------------
# 参数设置（可自行修改）
# ------------------------------------------------------------------------------------
m = 24  # 模态阶数
k = 31.0  # 无量纲波数
M = 0.5  # 马赫数

case = "a"  # ref[2]中figure1的case是a还是c
if case == "a":
    kmn = -44.2 + 1.1j  # Ref[2] a,b 轴向波数（可取复数）
else:
    kmn = -17.6 - 21.0j  # Ref[2] c,d 轴向波数（可取复数）


Z = 2.0 + 0.6j  # 壁面阻抗（复数）
omega = k  # 角频率，若使用不同无量纲化请自行修改
r_min = 1e-13  # 避免 r=0 的奇点，从一个很小的半径开始
r_max = 1.0
delta = 0.0002
isShearFlow = True  #
# ------------------------------------------------------------------------------------
# 将复数问题拆为实数 4 维系统
# y = [Re(p), Im(p), Re(dp/dr), Im(dp/dr)]
# ------------------------------------------------------------------------------------


# Ref[2] 公式(14)的对r的导数
def du0dr(r, delta):
    a = 1.0 - np.tanh(1.0 / delta)
    b = (1 + np.tanh(1.0 / delta)) / delta
    return M * (
        (1.0 - np.tanh((1.0 - r) / delta) ** 2) * (-1.0 / delta)
        + a * (-2.0 * (b + 1.0) * r + b)
    )


# Ref[1] 公式有误，使用Ref[2]的公式，应当*r，而不是+r
def u0(r, delta):
    return M * (
        np.tanh((1 - r) / delta)
        + (1 - r)
        * (1 - np.tanh(1 / delta))
        * ((1 + np.tanh(1 / delta)) / delta * r + (1 + r))
    )


def ode(r, y):
    p = y[0] + 1j * y[1]
    dp = y[2] + 1j * y[3]
    U = u0(r, delta) if isShearFlow else M
    coeff = (k - U * kmn) ** 2 - kmn**2

    if isShearFlow:
        denom = k - U * kmn
        dMdr = du0dr(r, delta)
        # 这里1e-12的作用？
        dMdrPart = np.where(np.abs(denom) > 1e-12, 2.0 * kmn / denom * dMdr, 0.0)
    else:
        dMdrPart = 0.0

    inv_r = np.where(r > 0.0, 1.0 / r, 0.0)
    inv_r2 = np.where(r > 0.0, m**2 / r**2, 0.0)
    d2p = -(inv_r + dMdrPart) * dp - (coeff - inv_r2) * p

    return np.vstack(
        [
            y[2],
            y[3],
            d2p.real,
            d2p.imag,
        ]
    )


def boundary_conditions(ya, yb):
    # r = r_min: 正则性 => r_min·p'(r_min) - m·p(r_min) = 0
    # 需要乘以r_min嘛？
    bc_center_real = r_min * ya[2] - m * ya[0]
    bc_center_imag = r_min * ya[3] - m * ya[1]

    # r = r_max: Ingard-Myers 边界条件
    p_wall = yb[0] + 1j * yb[1]
    dp_wall = yb[2] + 1j * yb[3]
    u_wall = u0(r_max, delta) if isShearFlow else M

    # 这个结合Ref[1] 公式(10),(19)组合与Ref[2]的公式(5)一致，都能得到如下边界条件
    impedance_condition = dp_wall - (k - u_wall * kmn) ** 2 / (1j * k * Z) * p_wall

    return np.array(
        [
            bc_center_real,
            bc_center_imag,
            impedance_condition.real,
            impedance_condition.imag,
        ]
    )


# ------------------------------------------------------------------------------------
# 初始猜测：使用 Bessel J_m(alpha r)，其中 alpha^2 = (k - M kmn)^2 - kmn^2
# ------------------------------------------------------------------------------------
alpha_squared = (k - M * kmn) ** 2 - kmn**2
alpha = np.sqrt(alpha_squared)

max_nodes = 1200000
r = np.linspace(r_min, r_max, max_nodes // 2)  # 这里点数太少导致幅值不对

p_guess = jv(m, alpha * r)
dp_guess = alpha * jvp(m, alpha * r)

y_guess = np.vstack(
    [
        p_guess.real,
        p_guess.imag,
        dp_guess.real,
        dp_guess.imag,
    ]
)


# ------------------------------------------------------------------------------------
# 调用 solve_bvp
# ------------------------------------------------------------------------------------
sol = solve_bvp(ode, boundary_conditions, r, y_guess, max_nodes=max_nodes)

if not sol.success:
    raise RuntimeError(f"solve_bvp 未收敛: {sol.message}")


# ------------------------------------------------------------------------------------
# 后处理与绘图
# ------------------------------------------------------------------------------------

# 这个根据Ref[2]要求归一化
p_wall = sol.y[0][-1] + 1j * sol.y[1][-1]
scale = 1.0 / p_wall
p_r = (sol.y[0] + 1j * sol.y[1]) * scale
dp_r = (sol.y[2] + 1j * sol.y[3]) * scale

print("p_wall (normalized to 1):", p_wall)

print("求解成功，最大残差 =", np.max(np.abs(sol.rms_residuals)))
# 读取文献数据
# 打开并读取文件
if case == "a":
    file_path = "./ref2_f1a.txt"  # Ref[2] a,b 数据文件路径
else:
    file_path = "./ref2_f1c.txt"  # Ref[2] c,d 数据文件路径
x_values = []  # 存储第一列数据
y_values = []  # 存储第二列数据

with open(file_path, "r") as file:
    for line in file:
        # 去掉行首行尾的空白字符，并按逗号分割
        x, y = line.strip().split(",")
        # 将字符串转换为浮点数，并添加到列表中
        x_values.append(float(x))
        y_values.append(float(y))
# 转换为numpy数组
x_values = np.array(x_values)
y_values = np.array(y_values)
print("x_values:", x_values)
print("y_values:", y_values)

plt.figure(figsize=(8, 5))
plt.plot(sol.x, np.abs(p_r), label="|p(r)|")
plt.plot(sol.x, p_r.real, "--", label="Re{p(r)}")
plt.plot(sol.x, p_r.imag, "--", label="Im{p(r)}")
plt.scatter(x_values, y_values, color="red", s=10, label="Ref[2] Fig.1 c data")
# p_r的幅值，相位模式
plt.xlabel("r")
plt.ylabel("声压")
plt.title(f"公式(13) BVP 解,m={m}")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
