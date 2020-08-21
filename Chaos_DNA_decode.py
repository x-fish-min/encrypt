import cv2 as cv
from scipy.integrate import odeint
from function import *
from pylab import *
mpl.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus']=False

lastKey = np.load('C:/Users/xym/Desktop/key.npy')
path = lastKey[0]
lastKey = lastKey[1:]
Key = list(map(float, lastKey))
I = cv.imread(path)
I1 = I[:, :, 0]
I2 = I[:, :, 1]
I3 = I[:, :, 2]
M = I.shape[0]
N = I.shape[1]
t = 4               #分块大小
SUM = M * N
u = 3.9999
xx0 = Key[0]
xx1 = Key[1]
ppx = np.zeros([M+1000])
ppy = np.zeros([N+1000])
# print(ppx)
ppx[0] = xx0
# print(ppx)
ppy[0] = xx1
for i in range(M+999):                      #进行M+999次循环，共得到M+1000点（包括初值）
    ppx[i+1] = u * ppx[i] * (1 - ppx[i])
for i in range(N+999):                      #进行M+999次循环，共得到N+1000点（包括初值）
    ppy[i+1] = u * ppy[i] * (1 - ppy[i])
ppx = ppx[1000 : len(ppx)]          #去除前1000点，获得更好的随机性
# print(ppx)
ppy = ppy[1000 : len(ppy)]

Ux = np.argsort(-ppx)
Uy = np.argsort(-ppy)

for i in range(N-1, -1, -1):
    temp = np.copy(I1[:, i])
    I1[:, i] = I1[:, Uy[i]]
    I1[:, Uy[i]] = temp
    temp = np.copy(I2[:, i])
    I2[:, i] = I2[:, Uy[i]]
    I2[:, Uy[i]] = temp
    temp = np.copy(I3[:, i])
    I3[:, i] = I3[:, Uy[i]]
    I3[:, Uy[i]] = temp

for i in range(M-1, -1, -1):
    temp = np.copy(I1[i, :])
    I1[i, :] = I1[Ux[i], :]
    I1[Ux[i], :] = temp
    temp = np.copy(I2[i, :])
    I2[i, :] = I2[Ux[i], :]
    I2[Ux[i], :] = temp
    temp = np.copy(I3[i, :])
    I3[i, :] = I3[Ux[i], :]
    I3[Ux[i], :] = temp

"""----------产生logistic混沌序列----------"""
x0 = Key[2]
p = np.zeros([SUM+1000])                                #预分配内存
p[0] = x0
for i in range(SUM+999):                                #进行SUM+999次循环，共得到SUM+1000点（包括初值）
    p[i+1] = u * p[i] * (1 - p[i])
p = p[1000 : len(p)]                                    #去除前1000点，获得更好的随机性

"""----------将p序列变换到0~255范围内整数，转换成M*N的二维矩阵R----------"""
p = mod(np.round(p * pow(10, 4)), 256)
R = p.reshape((M, N))                                   # 转成M行N列

"""----------求解混沌方程----------"""
# 求四个初值X0,Y0,Z0,H0
r = (M / t) * (N / t)
# X0=0.5008000000000001         #密钥敏感性测试
X0 = Key[3]
Y0 = Key[4]
Z0 = Key[5]
H0 = Key[6]
# X0=0.5056;                    #home图片
# Y0=0.505
# Z0=0.4564
# H0=0.3062
def func(w, time, m, n, p, q, r):
    # 给出位置矢量w，和三个参数m, n, p, q, r计算出
    # dx/dt, dy/dt, dz/dt的值
    x1, x2, x3, x4 = w
    # 直接与Chen超混沌系统的计算公式对应
    return np.array([m*(x2-x1), -x1*x3+q*x1+p*x2-x4, x1*x2-n*x3, x1+r])
# 创建时间点
time = np.arange(0, 500, 500/(r+3000))
# 调用odeint对Chen氏超混沌系统进行求解
A = odeint(func, (X0, Y0, Z0, H0), time, args=(36, 3, 28, -16, 0.42))

X = A[:, 0]
X = X[3000 : len(X)]        #去除前3000项，获得更好的随机性（求解陈氏系统的子函数多计算了3000点）
Y = A[:, 1]
Y = Y[3000 : len(Y)]
Z = A[:, 2]
Z = Z[3000 : len(Z)]
H = A[:, 3]
H = H[3000 : len(H)]

"""----------DNA编码----------"""
# X,Y分别决定I和R的DNA编码方式，有8种，1~8
X = mod(np.round(X * pow(10, 4)), 8) + 1
Y = mod(np.round(Y * pow(10, 4)), 8) + 1
Z = mod(np.round(Z * pow(10, 4)), 4)
Z[(Z == 0)] = 4                  # 加减法互换
Z[(Z == 1)] = 0
Z[(Z == 4)] = 1
H = mod(np.round(H * pow(10, 4)), 8) + 1
each_h = int(M / t)                     #每块的高边长
each_w = int(N / t)                     #每块的宽边长

Q_B = np.zeros([M, N])
Q_G = np.zeros([M, N])
Q_R = np.zeros([M, N])
for i in range(t*t, 1, -1):
    Q1_B = DNA_bian(fenkuai(t, I1, i), H[i-1])            # 对原始图像R通道每一个分块按X对应的序号进行DNA编码
    Q1_G = DNA_bian(fenkuai(t, I2, i), H[i-1])
    Q1_R = DNA_bian(fenkuai(t, I3, i), H[i-1])

    Q1_last_B = DNA_bian(fenkuai(t, I1, i-1), H[i-2])
    Q1_last_G = DNA_bian(fenkuai(t, I2, i-1), H[i-2])
    Q1_last_R = DNA_bian(fenkuai(t, I3, i-1), H[i-2])

    Q2_B = DNA_yunsuan(Q1_B, Q1_last_B, Z[i-1])           # 扩散前
    Q2_G = DNA_yunsuan(Q1_G, Q1_last_G, Z[i-1])
    Q2_R = DNA_yunsuan(Q1_R, Q1_last_R, Z[i-1])

    Q3 = DNA_bian(fenkuai(t, R, i), Y[i-1])             # 对R的每一个分块按Y对应的序号进行DNA编码

    Q4_B = DNA_yunsuan(Q2_B, Q3, Z[i-1])
    Q4_G = DNA_yunsuan(Q2_G, Q3, Z[i-1])
    Q4_R = DNA_yunsuan(Q2_R, Q3, Z[i-1])

    xx = int(floor(i / t) + 1)
    yy = int(mod(i, t))
    if yy == 0:
        xx = xx - 1
        yy = t
    Q_B[(xx - 1) * each_h : xx * each_h, (yy - 1) * each_w : yy * each_w] = DNA_jie(Q4_B, X[i-1])         # 将每一块合并成完整的图Q
    Q_G[(xx - 1) * each_h : xx * each_h, (yy - 1) * each_w : yy * each_w] = DNA_jie(Q4_G, X[i-1])
    Q_R[(xx - 1) * each_h : xx * each_h, (yy - 1) * each_w : yy * each_w] = DNA_jie(Q4_R, X[i-1])

Q5_B = DNA_bian(fenkuai(t, I1, 1), H[0])
Q5_G = DNA_bian(fenkuai(t, I2, 1), H[0])
Q5_R = DNA_bian(fenkuai(t, I3, 1), H[0])

Q6 = DNA_bian(fenkuai(t, R, 1), Y[0])

Q7_B = DNA_yunsuan(Q5_B, Q6, Z[0])
Q7_G = DNA_yunsuan(Q5_G, Q6, Z[0])
Q7_R = DNA_yunsuan(Q5_R, Q6, Z[0])

Q_B[0:each_h, 0:each_w] = DNA_jie(Q7_B, X[0])
Q_G[0:each_h, 0:each_w] = DNA_jie(Q7_G, X[0])
Q_R[0:each_h, 0:each_w] = DNA_jie(Q7_R, X[0])

Q_jiemi = np.zeros([M, N, 3], dtype=np.float32)
Q_jiemi[:, :, 0] = uint8(Q_R)
Q_jiemi[:, :, 1] = uint8(Q_G)
Q_jiemi[:, :, 2] = uint8(Q_B)

"""----------去除加密时补的零----------"""
t = Key[7]                # 加密时的分块参数，作为密钥
M1 = M / t
N1 = N / t
if M1 != 0 and N1 != 0:
    Q_jiemi = Q_jiemi[0 : int(M-t+M1), 0 : int(N-t+N1), :]
elif M1 != 0 and N1 == 0:
    Q_jiemi = Q_jiemi[0 : int(M-t+M1), 0 : N, :]
elif M1 == 0 and N1 != 0:
    Q_jiemi = Q_jiemi[0 : M, 0: int(N-t+N1), :]
else:
    Q_jiemi = Q_jiemi[0: M, 0: N, :]

cv.imwrite('./data/lena_d.png', cv.cvtColor(Q_jiemi, cv.COLOR_BGR2RGB))
plt.imshow(Q_jiemi.astype(np.uint8))
plt.show()
