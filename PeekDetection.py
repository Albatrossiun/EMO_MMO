'''''''''''''''''''''
# @FileName:PeekDetection.py
# @author:ZhaoXinYi
# @version:0.0.1
# @Date:2020.07.17
# @BSD
'''''''''''''''''''''
import numpy as np
import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
MAX = 9999999.0
MIN = -9999999.0

def distance(a, b):
    sum = 0 
    size = len(a)
    for i in range(size):
        sum += (a[i] - b[i])**2
    return sum**0.5

def ADP(D):
    # 构造距离矩阵
    size = len(D)
    dMatrix = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            if(i == j):
                dMatrix[i][j] = MAX
                continue
            dMatrix[j][i] = distance(D[i], D[j])
            dMatrix[i][j] = dMatrix[j][i]

    # 初始化临接距离数组
    neighboringDistance = []
    for i in range(size):
        tmp = MAX
        for j in range(size):
            if(dMatrix[i][j] < tmp):
                tmp = dMatrix[i][j]
        neighboringDistance.append(tmp)

    # P是返回的峰值检测结果
    P = []
    # 用来记录从D中删掉的点
    deletedPointSet = set(())
    # 当集合D不为空的时候
    while(size != len(deletedPointSet)):
        # 计算σ sigma 并更新临接距离数组
        sigma = MIN
        for i in range(size):
            tmpMin = MAX
            # 判断一下点i是不是已经被删除了
            if(i in deletedPointSet):
                continue
            for j in range(size):
                if(i == j):
                    continue
                # 判断一下点j是不是已经被删除了
                if(j in deletedPointSet):
                    continue
                if(dMatrix[i][j] < tmpMin):
                    tmpMin = dMatrix[i][j]
            neighboringDistance[i] = tmpMin
            if(tmpMin > sigma):
                sigma = tmpMin
        # print(sigma,"sigma")
        # 初始化峰值数组ψ Psi
        Psi = []
        
        for i in range(size):
            if (i in deletedPointSet):
                continue
            if (math.isclose(neighboringDistance[i],sigma)):
                Psi.append(i)
                deletedPointSet.add(i)
                break
        if(len(Psi) != 1):
            print("出现错误")
            return None
        # 遍历Psi数组
        i = 0
        while(True):
            PsiSize = len(Psi)
            if(i >= PsiSize):
                break
            index = Psi[i]
            for j in range(size):
                if(j == index):
                    continue
                # 如果j已经被删除了
                if(j in deletedPointSet):
                    continue
                # 如果这个点和 当前选中的Psi里的点i 之间的距离小于sigma 把这个点添加到Psi里边
                if(dMatrix[j][index] < sigma * 2):
                    Psi.append(j)
                    deletedPointSet.add(j)
            # 从Psi里拿出下一个点
            i = i + 1

        # 把下标转换成点
        tmp = []
        for i in range(len(Psi)):
            tmp.append(D[Psi[i]])
        # 把这一组峰值装进全部峰值集合里
        P.append(tmp.copy())
    return P

def PeekDection(PointList, valueList, nita = 0.03):
    D = []
    V = []
    max = MIN
    min = MAX
    for i in range(len(valueList)):
        if(valueList[i] > max):
            max = valueList[i]
        if(valueList[i] < min):
            min = valueList[i]
    limit = max - nita * (max - min)

    for i in range(len(valueList)):
        if(valueList[i] >= limit):
            D.append(PointList[i])
            V.append(valueList[i])
    allPeekSet = []
    while(len(D) != 0):
        print("峰值检测一次（二进制切割），检测的点的个数为[{}]".format(len(D)))
        peekSet = ADP(D)
        allPeekSet = allPeekSet + peekSet
        min = MAX
        ND = []
        NV = []
        for i in range(len(D)):
            if(V[i] < min):
                min = V[i]
        limit = (max + min) / 2
        for i in range(len(D)):
            if(V[i] > limit):
                ND.append(D[i])
                NV.append(V[i])
        D = ND
        V = NV
    return allPeekSet

def ThreeDImg(X,Y,Z):
    ax = plt.subplot(111, projection='3d')
    ax.scatter(X,Y,Z)
    plt.show()

def TwoDImg(X,Y):
    plt.scatter(X, Y, s=0.2)
    plt.show()

def testFunction(x):
    return math.cos(math.pi*0.3*x[0])**2 + math.sin(math.pi*0.3*x[1])**2

if __name__ == "__main__":
    X = []
    V = []
    for i in range(5000):
        a = 0 + random.random() * 9
        b = 0 + random.random() * 9
        X.append([a,b])
        V.append(testFunction(X[i]))
    ThreeDImg([x[0] for x in X],[x[1] for x in X],V)
    peekSet = PeekDection(X,V)
    TMP = []
    TV = []
    for i in range(len(peekSet)):
        for j in range(len(peekSet[i])):
            TMP.append(peekSet[i][j])
            TV.append(testFunction(peekSet[i][j]))
    ThreeDImg([x[0] for x in TMP], [x[1] for x in TMP], TV)
