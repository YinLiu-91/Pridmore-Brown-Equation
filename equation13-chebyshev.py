import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigs

# 定义Chebyshev微分矩阵
def chebyshev_diff_matrix(N):
    """生成N阶Chebyshev微分矩阵"""
    c = np.ones(N)
    c[0] = 2
    c[-1] = 2
    x = np.cos(np.pi * np.arange(N) / (N - 1))
    D = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i != j:
                D[i, j] = c[i] / c[j] / (x[i] - x[j])
            else:
                D[i, j] = -x[i] / (2 * (1 - x[i]**2))
    D = D.T  # 转置以匹配标准定义
    return D, x

# 定义Pridmore-Brown方程的系统矩阵
def pridmore_brown_matrix(N, k, m, M0, a):
    """生成Pridmore-Brown方程的系统矩阵"""
    D, x = chebyshev_diff_matrix(N)
    r = (1 + x) / 2  # 映射到[0, 1]
    
    # 马赫数分布
    M = M0 * (1 - r / a)
    dM_dr = -M0 / a
    
    # 系数矩阵
    A = 1 / r + (2 * k / (k - M * k)) * dM_dr
    B = (k - M * k)**2 - k**2 - m**2 / r**2
    
    # 系统矩阵
    L = D @ D + diags(A) @ D + diags(B)
    
    # 边界条件
    # 壁面（r=1）: \frac{d\hat{p}}{dr} + \alpha \hat{p} = 0
    # 中心（r=0）: \frac{d\hat{p}}{dr} = 0 (m=0) 或 \hat{p}(0) = 0 (m≠0)
    if m == 0:
        L[0, :] = D[0, :]  # 中心正则性条件
    else:
        L[0, :] = np.eye(N)[0, :]  # 中心正则性条件
    L[-1, :] = D[-1, :] + 100 * np.eye(N)[-1, :]  # 壁面阻抗条件（假设α=100）
    
    return L, r

# 参数设置
N = 50  # Chebyshev点数
k = 10  # 波数
m = 1   # azimuthal mode
M0 = 0.5  # 参考马赫数
a = 1    # 管道半径

# 生成系统矩阵
L, r = pridmore_brown_matrix(N, k, m, M0, a)

# 求解特征值问题
eigenvalues, eigenvectors = eigs(L, k=5, which='SM')  # 求最小的5个特征值

# 输出结果
print("特征值（轴向波数 k_{m,n}）:", eigenvalues)
print("特征向量（径向模态形状 \hat{p}(r)）:", eigenvectors)

# 可视化径向模态形状
import matplotlib.pyplot as plt

for i in range(eigenvectors.shape[1]):
    plt.plot(r, np.abs(eigenvectors[:, i]), label=f'Mode {i+1}')
plt.xlabel('r')
plt.ylabel('|p(r)|')
plt.legend()
plt.title('Radial Mode Shapes')
plt.show()