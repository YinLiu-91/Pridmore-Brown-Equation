import numpy as np
from scipy.integrate import solve_bvp
from scipy.special import jv, jvp
import matplotlib.pyplot as plt

# 参数
m = 24  # 模态阶数
k = 31.0  # 无量纲波数
M = 0.5  # 马赫数
kmn = 25.42 + 0.13j  # 已知轴向波数（来自前文计算）
Z = 3 - 0.5j  # 阻抗

# 定义一阶 ODE 系统
def ode(r, y):
    p, dpdr = y
    dp2dr2 = -dpdr / r - ((k - M * kmn)**2 - kmn**2 - m**2 / r**2) * p
    return np.vstack((dpdr, dp2dr2))

# 边界条件
def boundary_conditions(ya, yb):
    p0, dp0 = ya  # r=0 处的值
    p1, dp1 = yb  # r=1 处的值
    # r=0: 对称性 => dp/dr = 0
    bc1 = dp0
    # r=1: Ingard–Myers  =>  dp/dr = i(k-M*kmn)/Z * p
    bc2 = dp1 - 1j * (k - M * kmn) / Z * p1
    return np.array([bc1, bc2])

# 定义 r 范围
r = np.linspace(1e-6, 1, 200)  # 从 r=1e-6 到 r=1，共 200 个点

# 初始猜测
p_guess = jv(m, r)
dpdr_guess = jvp(m, r)
y_guess = np.vstack((p_guess, dpdr_guess))

# 使用 solve_bvp 求解边界值问题
sol = solve_bvp(ode, boundary_conditions, r, y_guess, tol=1e-5, max_nodes=1000)

# 检查求解结果
if sol.success:
    p_r = sol.y[0]
    dpdr_r = sol.y[1]
    print("求解成功")
else:
    print("求解失败：", sol.message)

# 绘图
plt.plot(r, np.abs(p_r), label=f'|p(r)|, m={m}')
plt.xlabel('r')
plt.ylabel('|p(r)|')
plt.title('声压沿径向 r 的分布（solve_bvp 方法）')
plt.legend()
plt.grid()
plt.show()