import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D
import math

MAX = 99999999.0
MIN = -99999999.0

colorList = ['k','r','y','g','c','b','m']

def distance(a,b,x,y):
    return ((a-x)**2 + (b-y)**2)**0.5

def distance1(a,b,x,y):
    return abs((a-x) + (b - y))

def ADP(D):
    # 构造距离矩阵
    size = len(D)
    dMatrix = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            if(i == j):
                dMatrix[i][j] = MAX
                continue
            dMatrix[j][i] = distance(D[i][0],D[i][1],D[j][0],D[j][1])
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
        print(sigma,"sigma")
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
                if(dMatrix[j][index] < sigma * 3):
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

#[x1,x2,y]
def PeekDection(PointList, nita = 0.1):
    D = []
    max = -1
    min = 9999
    for i in range(len(PointList)):
        if(PointList[i][2] > max):
            max = PointList[i][2]
        if(PointList[i][2] < min):
            min = PointList[i][2]
    limit = max - nita * (max - min)

    for i in range(len(PointList)):
        if(PointList[i][2] >= limit):
            D.append(PointList[i])
    print("从",len(D),"个点里")
    rrr = ADP(D)
    print("检测处",len(rrr),"个峰")
    ax = plt.subplot(111, projection='3d')
    for i in range(len(rrr)):
        X = []
        Y = []
        Z = []
        for j in range(len(rrr[i])):
            X.append(rrr[i][j][0])
            Y.append(rrr[i][j][1])
            Z.append(rrr[i][j][2])
        
        ax.scatter(X, Y, Z, c = colorList[i % len(colorList)])
    plt.show()


def ThreeDImg(X,Y,Z):
    ax = plt.subplot(111, projection='3d')
    ax.scatter(X,Y,Z)
    plt.show()

def TwoDImg(X,Y):
# 画散点图
    plt.scatter(X, Y, s=0.2)
    plt.show()



def objectionFunc(a,b):
    return math.cos(math.pi*0.3*a)**2 + math.sin(math.pi*0.3*b)**2


plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus'] = False
#matplotlib画图中中文显示会有问题，需要这两行设置默认字体
 
#plt.xlabel('X')
#plt.ylabel('Y')
#plt.xlim(xmax=9,xmin=0)
#plt.ylim(ymax=9,ymin=0)
#画两条（0-9）的坐标轴并设置轴标签x，y
 
x1 = np.random.normal(7,0.6,50) # 随机产生300个平均值为2，方差为1.2的浮点数
y1 = np.random.normal(7,0.6,50) # 随机产生300个平均值为2，方差为1.2的浮点数

x2 = np.random.normal(2,0.6,50)
y2 = np.random.normal(2,0.6,50)

x1 = np.hstack((x1,x2))
y1 = np.hstack((y1,y2))

x2 = np.random.normal(5,0.6,100)
y2 = np.random.normal(5,0.6,100)

x1 = np.hstack((x1,x2))
y1 = np.hstack((y1,y2))

x2 = np.random.normal(7,0.6,100)
y2 = np.random.normal(1,0.6,100)

x1 = np.hstack((x1,x2))
y1 = np.hstack((y1,y2))

X = []
Y = []
Z = []
pointList = []
for i in range(0):
    point = []
    point.append(x1[i])
    point.append(y1[i])
    point.append(objectionFunc(x1[i],y1[i]))
    X.append(point[0])
    Y.append(point[1])
    Z.append(point[2])
    pointList.append(point)
for i in range(5000):
    point = []
    a = 0 + random.random() * 9
    b = 0 + random.random() * 9
    point.append(a)
    point.append(b)
    point.append(objectionFunc(a,b))
    X.append(point[0])
    Y.append(point[1])
    Z.append(point[2])
    pointList.append(point)
'''
pointList = [[1,1,1],[1,1.2,1],[1,1.4,1],
[1.2,1,1],[1.2,1.2,1],[1.2,1.4,1],
[1.4,1,1],[1.4,1.2,1],[1.4,1.4,1],
[9,1,1],[9,1.2,1],[9,1.4,1],
[9.2,1,1],[9.2,1.2,1],[9.2,1.4,1],
[9.4,1,1],[9.4,1.2,1],[9.4,1.4,1]]
'''
PeekDection(pointList)

ThreeDImg(X,Y,Z)
TwoDImg(X,Y)


