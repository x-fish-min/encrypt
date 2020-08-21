import matplotlib.pyplot as plt
import cv2 as cv
from pylab import *
mpl.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus']=False

# 直方图分析
def Histogram(B, G, R):
    # 绘制原始图像各通道的直方图
    plt.hist(B.ravel(), 256)
    plt.title("原始图像B通道直方图")
    plt.xlabel("B通道像素灰度值")
    plt.ylabel("频数")
    plt.show()

    plt.hist(G.ravel(), 256)
    plt.title("原始图像G通道直方图")
    plt.xlabel("G通道像素灰度值")
    plt.ylabel("频数")
    plt.show()

    plt.hist(R.ravel(), 256)
    plt.title("原始图像R通道直方图")
    plt.xlabel("R通道像素灰度值")
    plt.ylabel("频数")
    plt.show()

# 相关性分析
def Correlation(B, G, R):
    # {
    # 先随机在0~H-1行和0~W-1列选中5000个像素点，
    # 计算水平相关性时，选择每个点的相邻的右边的点；
    # 计算垂直相关性时，选择每个点的相邻的下方的点；
    # 计算对角线相关性时，选择每个点的相邻的右下方的点。
    # }
    H, W = B.shape
    num = 5000       #随机取5000个像素点
    x = np.random.randint(0, H-1, 5000)      #生成5000个0~M-1的随机整数作为行
    y = np.random.randint(0, W-1, 5000)      #生成5000个0~N-1的随机整数作为列
    #预分配内存
    # 水平方向
    X_horizontal = np.zeros([num, 3])
    Y_horizontal = np.zeros([num, 3])
    # 垂直方向
    X_vertical = np.zeros([num, 3])
    Y_vertical = np.zeros([num, 3])
    # 对角线方向
    X_diagonal = np.zeros([num, 3])
    Y_diagonal = np.zeros([num, 3])

    for i in range(num):
        # 水平方向
        X_horizontal[i, 0] = B[x[i]][y[i]]
        Y_horizontal[i, 0] = B[x[i]+1][y[i]]
        X_horizontal[i, 1] = G[x[i]][y[i]]
        Y_horizontal[i, 1] = G[x[i]+1][y[i]]
        X_horizontal[i, 2] = R[x[i]][y[i]]
        Y_horizontal[i, 2] = R[x[i]+1][y[i]]
        # 垂直方向
        X_vertical[i, 0] = B[x[i]][y[i]]
        Y_vertical[i, 0] = B[x[i]][y[i]+1]
        X_vertical[i, 1] = G[x[i]][y[i]]
        Y_vertical[i, 1] = G[x[i]][y[i]+1]
        X_vertical[i, 2] = R[x[i]][y[i]]
        Y_vertical[i, 2] = R[x[i]][y[i]+1]
        # 对角线方向
        X_diagonal[i, 0] = B[x[i]][y[i]]
        Y_diagonal[i, 0] = B[x[i]+1][y[i]+1]
        X_diagonal[i, 1] = G[x[i]][y[i]]
        Y_diagonal[i, 1] = G[x[i]+1][y[i]+1]
        X_diagonal[i, 2] = R[x[i]][y[i]]
        Y_diagonal[i, 2] = R[x[i]+1][y[i]+1]

    """水平相关性绘图"""
    #B通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像B通道水平相邻元素相关性点图")
    plt.xlabel("B通道随机点像素灰度值")
    plt.ylabel("与该点相邻水平方向像素灰度值")
    #设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    #设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    #画散点图
    ax1.scatter(X_horizontal[:, 0], Y_horizontal[:, 0], c = 'b', marker = 'o')
    plt.show()

    # G通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像G通道水平相邻元素相关性点图")
    plt.xlabel("G通道随机点像素灰度值")
    plt.ylabel("与该点相邻水平方向像素灰度值")
    # 设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    # 设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    # 画散点图
    ax1.scatter(X_horizontal[:, 1], Y_horizontal[:, 1], c='b', marker='o')
    plt.show()

    # R通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像R通道水平相邻元素相关性点图")
    plt.xlabel("R通道随机点像素灰度值")
    plt.ylabel("与该点相邻水平方向像素灰度值")
    # 设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    # 设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    # 画散点图
    ax1.scatter(X_horizontal[:, 2], Y_horizontal[:, 2], c='b', marker='o')
    plt.show()

    """垂直相关性绘图"""
    # B通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像B通道垂直相邻元素相关性点图")
    plt.xlabel("B通道随机点像素灰度值")
    plt.ylabel("与该点相邻垂直方向像素灰度值")
    # 设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    # 设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    # 画散点图
    ax1.scatter(X_vertical[:, 0], Y_vertical[:, 0], c='b', marker='o')
    plt.show()

    # G通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像G通道垂直相邻元素相关性点图")
    plt.xlabel("G通道随机点像素灰度值")
    plt.ylabel("与该点相邻垂直方向像素灰度值")
    # 设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    # 设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    # 画散点图
    ax1.scatter(X_vertical[:, 1], Y_vertical[:, 1], c='b', marker='o')
    plt.show()

    # R通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像R通道垂直相邻元素相关性点图")
    plt.xlabel("R通道随机点像素灰度值")
    plt.ylabel("与该点相邻垂直方向像素灰度值")
    # 设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    # 设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    # 画散点图
    ax1.scatter(X_vertical[:, 2], Y_vertical[:, 2], c='b', marker='o')
    plt.show()

    """对角线相关性绘图"""
    # B通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像B通道对角线相邻元素相关性点图")
    plt.xlabel("B通道随机点像素灰度值")
    plt.ylabel("与该点相邻对角线方向像素灰度值")
    # 设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    # 设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    # 画散点图
    ax1.scatter(X_diagonal[:, 0], Y_diagonal[:, 0], c='b', marker='o')
    plt.show()

    # G通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像G通道对角线相邻元素相关性点图")
    plt.xlabel("G通道随机点像素灰度值")
    plt.ylabel("与该点相邻对角线方向像素灰度值")
    # 设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    # 设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    # 画散点图
    ax1.scatter(X_diagonal[:, 1], Y_diagonal[:, 1], c='b', marker='o')
    plt.show()

    # R通道
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("原始图像R通道对角线相邻元素相关性点图")
    plt.xlabel("R通道随机点像素灰度值")
    plt.ylabel("与该点相邻对角线方向像素灰度值")
    # 设置横、纵坐标范围
    plt.xlim(0, 255)
    plt.ylim(0, 255)
    # 设置横、纵坐标刻度值
    plt.xticks(range(0, 256, 15))
    plt.yticks(range(0, 256, 15))
    # 画散点图
    ax1.scatter(X_diagonal[:, 2], Y_diagonal[:, 2], c='b', marker='o')
    plt.show()

    """相关性计算"""
    # 水平相关性
    EX_B_h = np.mean(X_horizontal[:, 0])
    DX_B_h = np.std(X_horizontal[:, 0])
    EY_B_h = np.mean(Y_horizontal[:, 0])
    DY_B_h = np.std(Y_horizontal[:, 0])
    COV_B_h = np.sum((X_horizontal[:, 0] - EX_B_h) * (Y_horizontal[:, 0] - EY_B_h)) / num
    r_B_h = COV_B_h / (DX_B_h * DY_B_h)

    EX_G_h = np.mean(X_horizontal[:, 1])
    DX_G_h = np.std(X_horizontal[:, 1])
    EY_G_h = np.mean(Y_horizontal[:, 1])
    DY_G_h = np.std(Y_horizontal[:, 1])
    COV_G_h = np.sum((X_horizontal[:, 1] - EX_G_h) * (Y_horizontal[:, 1] - EY_G_h)) / num
    r_G_h = COV_G_h / (DX_G_h * DY_G_h)

    EX_R_h = np.mean(X_horizontal[:, 2])
    DX_R_h = np.std(X_horizontal[:, 2])
    EY_R_h = np.mean(Y_horizontal[:, 2])
    DY_R_h = np.std(Y_horizontal[:, 2])
    COV_R_h = np.sum((X_horizontal[:, 2] - EX_R_h) * (Y_horizontal[:, 2] - EY_R_h)) / num
    r_R_h = COV_R_h / (DX_R_h * DY_R_h)

    # 垂直相关性
    EX_B_v = np.mean(X_vertical[:, 0])
    DX_B_v = np.std(X_vertical[:, 0])
    EY_B_v = np.mean(Y_vertical[:, 0])
    DY_B_v = np.std(Y_vertical[:, 0])
    COV_B_v = np.sum((X_vertical[:, 0] - EX_B_v) * (Y_vertical[:, 0] - EY_B_v)) / num
    r_B_v = COV_B_v / (DX_B_v * DY_B_v)

    EX_G_v = np.mean(X_vertical[:, 1])
    DX_G_v = np.std(X_vertical[:, 1])
    EY_G_v = np.mean(Y_vertical[:, 1])
    DY_G_v = np.std(Y_vertical[:, 1])
    COV_G_v = np.sum((X_vertical[:, 1] - EX_G_v) * (Y_vertical[:, 1] - EY_G_v)) / num
    r_G_v = COV_G_v / (DX_G_v * DY_G_v)

    EX_R_v = np.mean(X_vertical[:, 2])
    DX_R_v = np.std(X_vertical[:, 2])
    EY_R_v = np.mean(Y_vertical[:, 2])
    DY_R_v = np.std(Y_vertical[:, 2])
    COV_R_v = np.sum((X_horizontal[:, 2] - EX_R_v) * (Y_horizontal[:, 2] - EY_R_v)) / num
    r_R_v = COV_R_v / (DX_R_v * DY_R_v)

    # 对角线相关性
    EX_B_d = np.mean(X_diagonal[:, 0])
    DX_B_d = np.std(X_diagonal[:, 0])
    EY_B_d = np.mean(Y_diagonal[:, 0])
    DY_B_d = np.std(Y_diagonal[:, 0])
    COV_B_d = np.sum((X_diagonal[:, 0] - EX_B_d) * (Y_diagonal[:, 0] - EY_B_d)) / num
    r_B_d = COV_B_d / (DX_B_d * DY_B_d)

    EX_G_d = np.mean(X_diagonal[:, 1])
    DX_G_d = np.std(X_diagonal[:, 1])
    EY_G_d = np.mean(Y_diagonal[:, 1])
    DY_G_d = np.std(Y_diagonal[:, 1])
    COV_G_d = np.sum((X_diagonal[:, 1] - EX_G_d) * (Y_diagonal[:, 1] - EY_G_d)) / num
    r_G_d = COV_G_d / (DX_G_d * DY_G_d)

    EX_R_d = np.mean(X_diagonal[:, 2])
    DX_R_d = np.std(X_diagonal[:, 2])
    EY_R_d = np.mean(Y_diagonal[:, 2])
    DY_R_d = np.std(Y_diagonal[:, 2])
    COV_R_d = np.sum((X_diagonal[:, 2] - EX_R_d) * (Y_horizontal[:, 2] - EY_R_d)) / num
    r_R_d = COV_R_d / (DX_R_d * DY_R_d)

    print('B通道水平相关系数：', r_B_h)
    print('G通道水平相关系数：', r_G_h)
    print('R通道水平相关系数：', r_R_h)
    print('B通道垂直相关系数：', r_B_v)
    print('G通道垂直相关系数：', r_G_v)
    print('R通道垂直相关系数：', r_R_v)
    print('B通道对角线相关系数：', r_B_d)
    print('G通道对角线相关系数：', r_G_d)
    print('R通道对角线相关系数：', r_R_d)

def Information_Entropy(B, G, R):
    H, W = B.shape
    SUM = H * W
    #B通道
    T_B = np.histogram(B.ravel(), 256, [0, 256])
    xxs_B=0
    #G通道
    T_G = np.histogram(G.ravel(), 256, [0, 256])
    xxs_G=0
    #R通道
    T_R = np.histogram(R.ravel(), 256, [0, 256])
    xxs_R=0

    for i in range(256):
        pp_B = T_B[0][i] / SUM                #每个灰度值占比，即每个灰度值的概率
        pp_G = T_G[0][i] / SUM
        pp_R = T_R[0][i] / SUM
        if pp_B != 0:
            xxs_B = xxs_B - pp_B * log2(pp_B)
        if pp_G != 0:
            xxs_G = xxs_G - pp_G * log2(pp_G)
        if pp_R != 0:
            xxs_R = xxs_R - pp_R * log2(pp_R)
    print('B通道信息熵：', xxs_B)
    print('G通道信息熵：', xxs_G)
    print('R通道信息熵：', xxs_R)

if __name__ == '__main__':
    img = cv.imread('./data/lena_e.png')
    B = img[:, :, 0]
    G = img[:, :, 1]
    R = img[:, :, 2]
    Histogram(B, G, R)
    Correlation(B, G, R)
    Information_Entropy(B, G, R)