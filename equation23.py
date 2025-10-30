
import numpy as np
from scipy.special import jv, jvp
from scipy.optimize import fsolve
from scipy.integrate import quad

# 参数
m = 24          # 模态阶数
k = 31.0        # 波数
M = 0.5         # 马赫数
Z = 3 - 0.5j    # 复阻抗
delta = 0.01    # 边界层厚度

# 流速剖面
def u0(r, delta):
    return M * (np.tanh((1 - r) / delta) + (1 - r) * (1 - np.tanh(1 / delta)) * ((1 + np.tanh(1 / delta)) / delta + r + (1 + r)))

# 计算 δI0 和 δI1
def integrand_i0(r, kmn, delta):
    return 1 - ((k - u0(r, delta) * kmn)**2) / ((k - M * kmn)**2)

def integrand_i1(r, kmn, delta):
    return 1 - ((k - M * kmn)**2) / ((k - u0(r, delta) * kmn)**2)

def delta_I0(kmn, delta):
    return quad(integrand_i0, 0, 1, args=(kmn, delta))[0]

def delta_I1(kmn, delta):
    return quad(integrand_i1, 0, 1, args=(kmn, delta))[0]

# 非线性方程
def equation(kmn):
    kmn = kmn[0] + 1j * kmn[1]
    alpha_squared = (k - M * kmn)**2 - kmn**2
    if alpha_squared.real < 0:
        return [1e6, 1e6]
    alpha = np.sqrt(alpha_squared)
    delta_i0 = delta_I0(kmn, delta)
    delta_i1 = delta_I1(kmn, delta)
    lhs = 1j * k * Z * (1 - (kmn**2 + m**2) * delta * delta_i1 * jv(m, alpha) / (alpha * jvp(m, alpha)))
    rhs = (k - M * kmn)**2 * (jvp(m, alpha) / (alpha * jv(m, alpha)) - delta * delta_i0)
    res = lhs - rhs
    return [res.real, res.imag]

# 初值
k0 = k / (1 + M)  # 均匀流近似
initial_guess = [k0, 0.1]

# 求解
sol = fsolve(equation, initial_guess)
kmn_solution = sol[0] + 1j * sol[1]
alpha_solution = np.sqrt((k - M * kmn_solution)**2 - kmn_solution**2)

print(f"公式(23)解得：")
print(f"k_mn = {kmn_solution:.4f}")
print(f"α = {alpha_solution:.4f}")
# 验证
alpha = alpha_solution
delta_i0 = delta_I0(kmn_solution, delta)
delta_i1 = delta_I1(kmn_solution, delta)
lhs = 1j * k * Z * (1 - (kmn_solution**2 + m**2) * delta * delta_i1 * jv(m, alpha) / (alpha * jvp(m, alpha)))
rhs = (k - M * kmn_solution)**2 * (jvp(m, alpha) / (alpha * jv(m, alpha)) - delta * delta_i0)
print("LHS:", lhs)
print("RHS:", rhs)
print("残差:", lhs - rhs)