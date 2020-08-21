import cv2 as cv
from scipy.integrate import odeint
from function import *
from pylab import *
mpl.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus']=False

path = './data/lena.jpg'    # 图片路径
I = cv.imread(path)
I1 = I[:, :, 0]             # B通道
I2 = I[:, :, 1]             # G通道
I3 = I[:, :, 2]             # R通道
t = 4                       # 分块大小(每行或每列t块)

"""---------加密过程----------"""
"""step1：补零"""
P1 = pad_0(I1, t)
P2 = pad_0(I2, t)
P3 = pad_0(I3, t)
M = P1.shape[0]
N = P1.shape[1]
SUM = M * N

"""step2：产生Logistic混沌序列"""
u=3.9999                                                #Logistic参数μ，自定为3.9999
x0 = (sum(P1[:]) + sum(P2[:])) / (255 * SUM * 2)        #计算得出Logistic初值x0
x0 = math.floor(x0 * pow(10, 4)) / pow(10, 4)           #保留4位小数
p = np.zeros([SUM+1000])                                #预分配内存
print("x0=", x0)
p[0] = x0
for i in range(SUM+999):                                #进行SUM+999次循环，共得到SUM+1000个元素值
    p[i+1] = u * p[i] * (1 - p[i])
p = p[1000 : len(p)]                                    #去除前1000个元素，获得更好的随机性

"""step3：将p序列变换到0~255范围内整数，转换成M*N的二维矩阵R"""
p = mod(np.round(p * pow(10, 4)), 256)
R = p.reshape((M, N))                                   #转成M行N列的随机矩阵R

"""step4：求解Chen氏超混沌系统"""
#求四个初值X0,Y0,Z0,H0
r = int((M / t) * (N / t))
X0 = sum(sum(P1 & 17)) / (17 * SUM)
Y0 = sum(sum(P2 & 34)) / (34 * SUM)
Z0 = sum(sum(P3 & 68)) / (68 * SUM)
H0 = sum(sum(P1 & 136)) / (136 * SUM)
#保留四位小数
X0 = round(X0, 4)
Y0 = round(Y0, 4)
Z0 = round(Z0, 4)
H0 = round(H0, 4)
print("X0=", X0)
print("Y0=", Y0)
print("Z0=", Z0)
print("H0=", H0)
#根据初值，求解Chen氏超混沌系统，得到四个混沌序列
def func(w, time, m, n, p, q, r):
    # 给出位置矢量w，和三个参数m, n, p, q, r计算出
    # dx/dt, dy/dt, dz/dt的值
    x1, x2, x3, x4 = w
    # 直接与Chen超混沌系统的计算公式对应
    return np.array([m*(x2-x1), -x1*x3+q*x1+p*x2-x4, x1*x2-n*x3, x1+r])
# 创建时间点
time = np.arange(0, 500, 500 / (r + 3000))
# 调用odeint对Chen氏超混沌系统进行求解
A = odeint(func, (X0, Y0, Z0, H0), time, args=(36, 3, 28, -16, 0.42))

X = A[:, 0]
X = X[3000 : len(X)]        # 去除前3000项，获得更好的随机性（求解陈氏系统的子函数多计算了3000点）
Y = A[:, 1]
Y = Y[3000 : len(Y)]
Z = A[:, 2]
Z = Z[3000 : len(Z)]
H = A[:, 3]
H = H[3000 : len(H)]

"""step5：DNA编码"""
# X,Y分别决定P和R的DNA编码方式，有8种，1~8
# Z决定运算方式，有4种，0~3，0表示加，1表示减，2表示异或，3表示同或
# H表示DNA解码方式，有8种，1~8
X = mod(np.round(X * pow(10, 4)), 8) + 1
Y = mod(np.round(Y * pow(10, 4)), 8) + 1
Z = mod(np.round(Z * pow(10, 4)), 4)
H = mod(np.round(H * pow(10, 4)), 8) + 1
each_h = int(M / t)                     #每块的高边长
each_w = int(N / t)                     #每块的宽边长

Q_B = np.zeros([M, N])
Q_G = np.zeros([M, N])
Q_R = np.zeros([M, N])

Q2 = DNA_bian(fenkuai(t, R, 1), Y[0])

# B通道
Q1_B = DNA_bian(fenkuai(t, P1, 1), X[0])
Q_last_B = DNA_yunsuan(Q1_B, Q2, Z[0])
Q_B[0:each_h, 0:each_w] = DNA_jie(Q_last_B, H[0])
# G通道
Q1_G = DNA_bian(fenkuai(t, P2, 1), X[0])
Q_last_G = DNA_yunsuan(Q1_G, Q2, Z[0])
Q_G[0:each_h, 0:each_w] = DNA_jie(Q_last_G, H[0])
# R通道
Q1_R = DNA_bian(fenkuai(t, P3, 1), X[0])
Q_last_R = DNA_yunsuan(Q1_R, Q2, Z[0])
Q_R[0:each_h, 0:each_w] = DNA_jie(Q_last_R, H[0])

for i in range(2, t*t+1):
    Q1_B = DNA_bian(fenkuai(t, P1, i), X[i-1])        # 对原始图像R通道每一个分块按X对应的序号进行DNA编码
    Q1_G = DNA_bian(fenkuai(t, P2, i), X[i-1])
    Q1_R = DNA_bian(fenkuai(t, P3, i), X[i-1])

    Q2 = DNA_bian(fenkuai(t, R, i), Y[i-1])           # 对R的每一个分块按Y对应的序号进行DNA编码
    # B通道
    Q3_B = DNA_yunsuan(Q1_B, Q2, Z[i-1])
    Q4_B = DNA_yunsuan(Q3_B, Q_last_B, Z[i-1])
    Q_last_B = Q4_B
    # G通道
    Q3_G = DNA_yunsuan(Q1_G, Q2, Z[i-1])
    Q4_G = DNA_yunsuan(Q3_G, Q_last_G, Z[i-1])
    Q_last_G = Q4_G
    # R通道
    Q3_R = DNA_yunsuan(Q1_R, Q2, Z[i-1])              # 对上面两个编码好的块按Z对应的序号进行DNA运算
    Q4_R = DNA_yunsuan(Q3_R, Q_last_R, Z[i-1])        # 运算结果在和前一块按Z对应的序号再一次进行运算，称为扩散
    Q_last_R = Q4_R

    xx = int(floor(i / t) + 1)
    yy = int(mod(i, t))
    if yy == 0:
        xx = xx - 1
        yy = t
    Q_B[(xx - 1) * each_h : xx * each_h, (yy - 1) * each_w : yy * each_w] = DNA_jie(Q4_B, H[i-1])         # 将每一块合并成完整的图Q
    Q_G[(xx - 1) * each_h : xx * each_h, (yy - 1) * each_w : yy * each_w] = DNA_jie(Q4_G , H[i-1])
    Q_R[(xx - 1) * each_h : xx * each_h, (yy - 1) * each_w : yy * each_w] = DNA_jie(Q4_R, H[i-1])

Q_B = uint8(Q_B)
Q_G = uint8(Q_G)
Q_R = uint8(Q_R)

"""----------抗裁剪----------"""
xx0 = sum(P2[:]) / (255 * SUM)                      # G通道：平均灰度值，作为密钥
xx0 = math.floor(xx0 * pow(10, 4)) / pow(10, 4)          # 保留4位小数
xx1 = sum(P1[:]) / (255 * SUM)                      # B通道：平均灰度值，作为密钥
xx1 = math.floor(xx1 * pow(10, 4)) / pow(10, 4)          # 保留4位小数
ppx = zeros([M+1000])        # 预分配内存
ppy = zeros([N+1000])
print("xx0=", xx0)
print("xx1=", xx1)
ppx[0] = xx0
ppy[0] = xx1
for i in range(M+999):                      # 进行SUM+999次循环，共得到SUM+1000点（包括初值）
    ppx[i+1] = u * ppx[i] * (1-ppx[i])
for i in range(N+999):                      # 进行SUM+999次循环，共得到SUM+1000点（包括初值）
    ppy[i+1] = u * ppy[i] * (1-ppy[i])
ppx = ppx[1000 : len(ppx)]                  # 去除前1000点，获得更好的随机性
ppy = ppy[1000 : len(ppy)]

Ux = np.argsort(-ppx)
Uy = np.argsort(-ppy)

for i in range(M):
    temp = np.copy(Q_B[i, :])
    Q_B[i, :] = Q_B[Ux[i], :]
    Q_B[Ux[i], :] = temp
    temp = np.copy(Q_G[i, :])
    Q_G[i, :] = Q_G[Ux[i], :]
    Q_G[Ux[i], :] = temp
    temp = np.copy(Q_R[i, :])
    Q_R[i, :] = Q_R[Ux[i], :]
    Q_R[Ux[i], :] = temp

for i in range(N):
    temp = np.copy(Q_B[:, i])
    Q_B[:, i] = Q_B[:, Uy[i]]
    Q_B[:, Uy[i]] = temp
    temp = np.copy(Q_G[:, i])
    Q_G[:, i] = Q_G[:, Uy[i]]
    Q_G[:, Uy[i]] = temp
    temp = np.copy(Q_R[:, i])
    Q_R[:, i] = Q_R[:, Uy[i]]
    Q_R[:, Uy[i]] = temp

Q_jiami = np.zeros([M, N, 3], dtype=np.float32)
Q_jiami[:, :, 0] = Q_B
Q_jiami[:, :, 1] = Q_G
Q_jiami[:, :, 2] = Q_R

savePath = './data/lena_e.png'
cv.imwrite(savePath, Q_jiami)
plt.imshow(Q_jiami.astype(np.uint8))
plt.show()

theKey = [savePath, xx0, xx1, x0, X0, Y0, Z0, H0, t]
np.save('./data/key', theKey)

