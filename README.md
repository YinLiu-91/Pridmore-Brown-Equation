# 说明
1. equation13-solve-bvp.py此文件根据给定的k，kmn等生成沿r的p函数，并与文献对比
1. equation22.py求解带流动的Pridmore-Brown方程的轴向波数kmn，考虑了brambley边界
1. 可以用equation22.py求kmn，将kmn带入equation13-solve-bvp.py中，求解p(r)分部
1. equation22.py中，设置如下初始值，可以得到ref[22]中的figure1的解：
   - initial_guess = [-44.2 ,1.1] # ref[2] figure 1a
   - initial_guess = [-17.6 ,- 21.0] # ref[2] figure 1c
1. 要想使用equation22.py给出某个(m,n)下正确的解，给出初始值是关键，不然只给m，对应的n是不确定的，可以利用理论方法如liudaren给出的代码中求出初始值，再使用这个脚本求解
# 参考文献
Ref:[1] On the prediction of far-field fan noise attenuation due to liners considering uniform and shear flows
Ref:[2] A well-posed boundary condition for acoustic liners in straight ducts with flow
Ref:[3] Impedance boundary condition in frequency domain DG code  for modelling liner effects in complicated intake ducts
