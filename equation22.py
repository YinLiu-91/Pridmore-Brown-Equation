# 文献 ：On the prediction of far-field fan noise attenuation due to liners
#       considering uniform and shear flows
import numpy as np
from scipy.special import jv, jvp
from scipy.optimize import fsolve
from scipy.integrate import quad

# 参数
m = 24          # 模态阶数
k = 31.0        # 波数
# M = -0.5         # 马赫数 Ref[3]
M = 0.5         # 马赫数 Ref[1][2]
# Z = 3 - 0.5j    # 复阻抗 Ref[3]
Z = 2 + 0.6j    # 复阻抗 Ref[1][2]
# delta = 0.02    # 边界层厚度
delta = 0.0002    # 边界层厚度
a=1.0          # 管道半径
# 流速剖面,Ref[1],Ref[3]都是这个形式
def u0(r, delta):
    return M * (np.tanh((1 - r) / delta) + (1 - r) * (1 - np.tanh(1 / delta)) * ((1 + np.tanh(1 / delta)) / delta + r + (1 + r)))
#ref[2]中的速度剖面
def u01(r, delta):
    return M * (
        np.tanh((1 - r) / delta)
        + (1 - r)
        * (1 - np.tanh(1 / delta))
        * ((1 + np.tanh(1 / delta)) / delta * r + (1 + r))
    )
# 计算 δI0 和 δI1
def integrand_i0(r, kmn, delta):
    return 1 - ((k - u0(r, delta) * kmn)**2) / ((k - M * kmn)**2)

def integrand_i1(r, kmn, delta):
    return 1 - ((k - M * kmn)**2) / ((k - u0(r, delta) * kmn)**2)

def delta_I0(kmn, delta):
    return 0.0 if delta==0.0 else  quad(integrand_i0, 0, a, args=(kmn, delta))[0]

def delta_I1(kmn, delta):
    return 0.0 if delta==0.0 else quad(integrand_i1, 0, a, args=(kmn, delta))[0]

# 非线性方程
def equation(kmn):
    kmn = kmn[0] + 1j * kmn[1]
    alpha_squared = (k - M * kmn)**2 - kmn**2
    if alpha_squared.real < 0:
        pass
        print("alpha_squared.real < 0")
        # return [1e6, 1e6]
    alpha = np.sqrt(alpha_squared)
    delta_i0 = delta_I0(kmn, delta)
    delta_i1 = delta_I1(kmn, delta)
    lhs = 1j * k * Z * (1 - (kmn**2 + m**2) * delta * delta_i1 * jv(m, alpha) / (alpha * jvp(m, alpha)))
    # 公式22中除数中有kmn,但是交叉验证后为alpha,(见A well-posed boundary condition for acoustic liners in
    # straight ducts with flow文中公式8下面的公式，同除以alphaJ‘),并且下面为jvp，上面为jv
    print("delta_i0:", delta_i0)
    print("delta_i1:", delta_i1)
    rhs = (k - M * kmn)**2 * (jv(m, alpha) / (alpha * jvp(m, alpha)) - delta * delta_i0)
    res = lhs - rhs
    return [res.real, res.imag]

# 初值
k0 = k / (1 + M)  # 均匀流近似
# initial_guess = [k0, 0.1]
# initial_guess = [7.2278,-2.4025]

# 下面两个ref[2]中的kmn都作为初值，都能得到与初始值相近的解
initial_guess = [-44.2 , 1.1] # ref[2] figure 1a
# initial_guess = [-20 ,-0.1] # ref[2] figure 1a
# initial_guess = [-17.6 ,- 21.0] # ref[2] figure 1c

# 求解
sol = fsolve(equation, initial_guess)
print("sol:", sol)
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
# 
rhs = (k - M * kmn_solution)**2 * (jv(m, alpha) / (alpha * jvp(m, alpha)) - delta * delta_i0)
print("LHS:", lhs)
print("RHS:", rhs)
print("残差:", lhs - rhs)

Amn=1.0

# 给出沿着r与x方向上的压力分部
r=np.linspace(0,a,200)
x=np.linspace(0,1,100)
# z
p=Amn*jv(m,alpha_solution*(a-r[:,np.newaxis]))*np.exp(1j*kmn_solution*x[np.newaxis,:])
import matplotlib.pyplot as plt
plt.imshow(np.real(p),extent=[0,1,0,a],aspect='auto')
plt.colorbar()
plt.title('Pressure Distribution Re(p) along r and x')
plt.xlabel('x')
plt.ylabel('r')
plt.show()

# 取r=1处沿x方向的压力分布
x_line = np.linspace(0, 1, 100)
print("p shape:", p.shape)
p_line=p[0,:]
# p_line = Amn * jvp(m, alpha_solution * (a )) * np.exp(1j * kmn_solution * x_line)
plt.plot(x_line, np.abs(p_line))
plt.title('Pressure Distribution Re(p) at r=1 along x')
plt.xlabel('x')
plt.ylabel('Re(p)')
plt.show()