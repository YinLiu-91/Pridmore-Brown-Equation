import numpy as np
from scipy.special import jv, jvp
from scipy.optimize import fsolve

# 参数
m = 24          # 模态阶数
k = 31.0        # 波数
M = 0.5         # 马赫数
Z = 3 - 0.5j    # 复阻抗

# 非线性方程（公式20）
def equation(kmn):
    kmn = kmn[0] + 1j*kmn[1]
    alpha_squared = (k - M*kmn)**2 - kmn**2
    if alpha_squared.real < 0:
        return [1e6, 1e6]
    alpha = np.sqrt(alpha_squared)
    lhs = alpha * jvp(m, alpha) / jv(m, alpha)
    rhs = -1j * k / Z * (1 - M * kmn / k)**2
    res = lhs - rhs
    return [res.real, res.imag]

# 初值
k0 = k / (1 + M)  # 均匀流近似
initial_guess = [k0, 0.1]

# 求解
sol = fsolve(equation, initial_guess)
kmn_solution = sol[0] + 1j*sol[1]
alpha_solution = np.sqrt((k - M*kmn_solution)**2 - kmn_solution**2)

print(f"公式(20)解得：")
print(f"k_mn = {kmn_solution:.4f}")
print(f"α = {alpha_solution:.4f}")

# 验证
alpha = alpha_solution
lhs = alpha * jvp(m, alpha) / jv(m, alpha)
rhs = -1j * k / Z * (1 - M * kmn_solution / k)**2
print("LHS:", lhs)
print("RHS:", rhs)
print("残差:", lhs - rhs)

