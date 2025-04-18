import numpy, xarray, pickle
from numpy import linalg as LA
from LIM_utils import LIM, Nyquist_check, Q_test, heatmap, Error_test, tau_test
import matplotlib.pyplot as plt

tau0 = 2
ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"

x = xarray.open_dataarray(ddir+"/state_vectors.nc")

#因为平滑的原因，前两天和后两天没有数据，所以舍去。
x = x.to_numpy()
print(x.shape)

output = LIM(x,tau0)
L = output['L']
Q = output['Q']

decayT = output['decayT']
frequency = output['frequency']

#正交性测试 如果正交则应该为0
#norm1 = numpy.linalg.norm(L.T@L - L@L.T, ord=2)
#norm2 = numpy.linalg.norm(L, ord=2)
#print(norm1/(norm2*norm2))

#Q_test: 矩阵Q的特征值应为正
#Q_test(Q,Q_plot='no')
#Nyquest_test: 避免出现频谱泄露或者混淆现象
#Nyquist_check(x,tau0)
#tau_test: 确保矩阵L同tau的选择无关
#tau_test(x,tau0,50)

Error_test(x,tau0,unit_days=121,split_num=30)

#保存LIM 模型的参数
with open(ddir+'/LIM_params.pkl', 'wb') as f:
    pickle.dump(output, f)