from pylab import *

# 补零函数
def pad_0(img, t):
    M = img.shape[0]
    N = img.shape[1]
    M1 = int(ceil(M / t))
    N1 = int(ceil(N / t))
    P = np.zeros([M1 * t, N1 * t], dtype=np.uint8)
    P[0 : M, 0 : N] = img[:, :]
    return P


# 分块函数，t为分块数量，I为要分块的图像，num为返回第几大块
def fenkuai(t, array, num):
    h = array.shape[0]
    w = array.shape[1]
    eh = h / t
    ew = w / t
    x = int(math.floor(num / t) + 1)         #第几大行
    y = int(np.mod(num, t))                  #第几大列
    if y == 0:
        x = x - 1
        y = t
    f = array[int(eh*(x-1)) : int(eh*x), int(ew*(y-1)) : int(ew*y)]
    return f


# 子函数 DNA编码
def DNA_bian(array, num):
    h_axi = array.shape[0]
    w_axi = array.shape[1]
    a1 = (array.astype(np.int8) & 192) / 64
    a2 = (array.astype(np.int8) & 48) / 16
    a3 = (array.astype(np.int8) & 12) / 4
    a4 = (array.astype(np.int8) & 3)
    A = np.hstack((a1, a2, a3, a4))
    fv = np.zeros([h_axi, 4*w_axi], dtype=np.str_)
    if num == 1 or num == 2:
        for i in range(h_axi):
            for j in range(4*w_axi):
                if A[i][j] == 0:
                    fv[i, j] = 'A'
                elif A[i][j] == 3:
                    fv[i, j] = 'T'
                elif A[i][j] == 2:
                    if num ==1:
                        fv[i, j] = 'C'
                    else:
                        fv[i, j] = 'G'
                else:
                    if num == 1:
                        fv[i, j] = 'G'
                    else:
                        fv[i, j] = 'C'
    elif num == 3 or num == 4:
        for i in range(h_axi):
            for j in range(4*w_axi):
                if A[i][j] == 1:
                    fv[i, j] = 'A'
                elif A[i][j] == 2:
                    fv[i, j] = 'T'
                elif A[i][j] == 0:
                    if num == 3:
                        fv[i, j] = 'G'
                    else:
                        fv[i, j] = 'C'
                else:
                    if num == 3:
                        fv[i, j] = 'C'
                    else:
                        fv[i, j] = 'G'
    elif num == 5 or num == 6:
        for i in range(h_axi):
            for j in range(4*w_axi):
                if A[i][j] == 2:
                    fv[i, j] = 'A'
                elif A[i][j] == 1:
                    fv[i, j] = 'T'
                elif A[i][j] == 0:
                    if num == 5:
                        fv[i, j] = 'G'
                    else:
                        fv[i, j] = 'C'
                else:
                    if num == 5:
                        fv[i, j] = 'C'
                    else:
                        fv[i, j] = 'G'
    else:
        for i in range(h_axi):
            for j in range(4*w_axi):
                if A[i][j] == 3:
                    fv[i, j] = 'A'
                elif A[i][j] == 0:
                    fv[i, j] = 'T'
                elif A[i][j] == 1:
                    if num == 7:
                        fv[i, j] = 'G'
                    else:
                        fv[i, j] = 'C'
                else:
                    if num == 7:
                        fv[i, j] = 'C'
                    else:
                        fv[i, j] = 'G'
    return fv


# DNA运算
# num为0时，表示加；num为1时，表示减；num为2时，表示异或；num为3时表示同或
def DNA_yunsuan(arr1, arr2, num):
    m = arr1.shape[0]
    n = arr1.shape[1]
    fv = np.zeros([m, n], dtype=np.str_)            # 预分配内存
    if num == 0:                     # 加法
        for i in range(m):
            for j in range(n):
                if arr1[i, j] == 'A':
                    fv[i, j] = arr2[i, j]
                elif arr1[i, j] == 'T':
                    if arr2[i, j] == 'A':
                        fv[i, j] = 'T'
                    elif arr2[i, j] == 'T':
                        fv[i, j] = 'C'
                    elif arr2[i, j] == 'C':
                        fv[i, j] = 'G'
                    else:
                        fv[i, j] = 'A'
                elif arr1[i, j] == 'C':
                    if arr2[i, j] == 'A':
                        fv[i, j] = 'C'
                    elif arr2[i, j] == 'T':
                        fv[i, j] = 'G'
                    elif arr2[i, j] == 'C':
                        fv[i, j] = 'A'
                    else:
                        fv[i, j] = 'T'
                else:
                    if arr2[i, j] == 'A':
                        fv[i, j] = 'G'
                    elif arr2[i, j] == 'T':
                        fv[i, j] = 'A'
                    elif arr2[i, j] == 'C':
                        fv[i, j] = 'T'
                    else:
                        fv[i, j] = 'C'

    elif num == 1:                   # 减法
        for i in range(m):
            for j in range(n):
                if arr2[i, j] == 'A':
                    fv[i, j] = arr1[i, j]
                elif arr2[i, j] == 'T':
                    if arr1[i, j] == 'A':
                        fv[i, j] = 'G'
                    elif arr1[i, j] == 'T':
                        fv[i, j] = 'A'
                    elif arr1[i, j] == 'C':
                        fv[i, j] = 'T'
                    else:
                        fv[i, j] = 'C'
                elif arr2[i, j] == 'C':
                    if arr1[i, j] == 'A':
                        fv[i, j] = 'C'
                    elif arr1[i, j] == 'T':
                        fv[i, j] = 'G'
                    elif arr1[i, j] == 'C':
                        fv[i, j] = 'A'
                    else:
                        fv[i, j] = 'T'
                else:
                    if arr1[i, j] == 'A':
                        fv[i, j] = 'T'
                    elif arr1[i, j] == 'T':
                        fv[i, j] = 'C'
                    elif arr1[i, j] == 'C':
                        fv[i, j] = 'G'
                    else:
                        fv[i, j] = 'A'

    elif num == 2:                   # 异或
        for i in range(m):
            for j in range(n):
                if arr1[i, j] == arr2[i, j]:
                    fv[i, j] = 'C'
                elif (arr1[i, j] == 'T' and arr2[i, j] == 'A') or (arr1[i, j] == 'A' and arr2[i, j] == 'T') or (
                      arr1[i, j] == 'G' and arr2[i, j] == 'C') or (arr1[i, j] == 'C' and arr2[i, j] == 'G'):
                    fv[i, j] = 'G'
                elif (arr1[i, j] == 'C' and arr2[i, j] == 'A') or (arr1[i, j] == 'A' and arr2[i, j] == 'C') or (
                      arr1[i, j] == 'G' and arr2[i, j] == 'T') or (arr1[i, j] == 'T' and arr2[i, j] == 'G'):
                    fv[i, j] = 'A'
                else:
                    fv[i, j] = 'T'

    else:                            # 同或
        for i in range(m):
            for j in range(n):
                if arr1[i, j] == arr2[i, j]:
                    fv[i, j] = 'A'
                elif (arr1[i, j] == 'T' and arr2[i, j] == 'A') or (arr1[i, j] == 'A' and arr2[i, j] == 'T') or (
                      arr1[i, j] == 'G' and arr2[i, j] == 'C') or (arr1[i, j] == 'C' and arr2[i, j] == 'G'):
                    fv[i, j] = 'T'
                elif (arr1[i, j] == 'C' and arr2[i, j] == 'A') or (arr1[i, j] == 'A' and arr2[i, j] == 'C') or (
                      arr1[i, j] == 'G' and arr2[i, j] == 'T') or (arr1[i, j] == 'T' and arr2[i, j] == 'G'):
                    fv[i, j] = 'C'
                else:
                    fv[i, j] = 'G'
    return fv


# DNA解码
def DNA_jie(array, num):
    m = array.shape[0]
    n = array.shape[1]
    e = int(n / 4)
    Out = np.zeros([m, n])
    if num ==1 or num == 2:
        for i in range(m):
            for j in range(n):
                if array[i, j] == 'A':
                    Out[i, j] = 0
                elif array[i, j] == 'T':
                    Out[i, j] = 3
                elif array[i, j] == 'G':
                    if num == 1:
                        Out[i, j] = 1
                    else:
                        Out[i, j] = 2
                else:
                    if num == 1:
                        Out[i, j] = 2
                    else:
                        Out[i, j] = 1
    elif num == 3 or num == 4:
        for i in range(m):
            for j in range(n):
                if array[i, j] == 'A':
                    Out[i, j] = 1
                elif array[i, j] == 'T':
                    Out[i, j] = 2
                elif array[i, j] == 'G':
                    if num == 3:
                        Out[i, j] = 0
                    else:
                        Out[i, j] = 3
                else:
                    if num == 3:
                        Out[i, j] = 3
                    else:
                        Out[i, j] = 0

    elif num == 5 or num == 6:
        for i in range(m):
            for j in range(n):
                if array[i, j] == 'A':
                    Out[i, j] = 2
                elif array[i, j] == 'T':
                    Out[i, j] = 1
                elif array[i, j] == 'G':
                    if num == 5:
                        Out[i, j] = 0
                    else:
                        Out[i, j] = 3
                else:
                    if num == 5:
                        Out[i, j] = 3
                    else:
                        Out[i, j] = 0

    else:               # num == 7 | | num == 8
        for i in range(m):
            for j in range(n):
                if array[i, j] == 'A':
                    Out[i, j] = 3
                elif array[i, j] == 'T':
                    Out[i, j] = 0
                elif array[i, j] == 'G':
                    if num == 7:
                        Out[i, j] = 1
                    else:
                        Out[i, j] = 2
                else:
                    if num == 7:
                        Out[i, j] = 2
                    else:
                        Out[i, j] = 1

    O1 = Out[:, 0: e]
    O2 = Out[:, e : 2 * e]
    O3 = Out[:, 2 * e : 3 * e]
    O4 = Out[:, 3 * e : 4 * e]
    fv = O1 * 64 + O2 * 16 + O3 * 4 + O4
    return fv

