# coding=utf-8
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## 单个个体
# 对于单个个体而言，每个个体有多个维度，每个维度有一定的取值区间，对于精度也有一定的要求
# 另外 每个个体应当保存自己的适应度

MAX = 99999999

class Gene:
    def __init__(self, limit, precision):
        self._limit = limit
        self._precision = precision
        # 计算二进制串的位数
        # [1,100] 精度为10
        space = (limit[1] - limit[0]) * pow(10.0, precision) # 自变量不同取值的个数
        len = 1
        while(pow(2.0, len) < space):
            len = len + 1
        self._len = len
        self._binaryGrayList = []

    def Copy(self):
        newGene = Gene(self._limit, self._precision)
        newGene._len = self._len
        newGene._binaryGrayList = []
        for i in range(self._len):
            newGene._binaryGrayList.append(self._binaryGrayList[i])
        return newGene       

    def Random(self):
        self._binaryGrayList = []
        for _ in range(self._len):
            tmp = 1
            if random.random() > 0.5:
                tmp = 0
            self._binaryGrayList.append(tmp)

    def Decode(self):
        binaryList = []
        tmp = self._binaryGrayList[0]
        binaryList.append(tmp)
        for i in range(self._len - 1):
            tmp = self._binaryGrayList[i+1]
            tmp1 = binaryList[i] ^ tmp
            binaryList.append(tmp1)
        sum = 0
        for i in range(self._len):
            sum = sum * 2 + binaryList[i]
        sum = sum * 1.0 / (pow(2, self._len) - 1)
        return self._limit[0] + sum * (self._limit[1] - self._limit[0]) 

    def GetLen(self):
        return self._len

    def CrossOver(self, gene, index):
        '''
        tmp = gene._binaryGrayList[index:]
        i = index
        while(i < self._len):
            gene._binaryGrayList[i] = self._binaryGrayList[i]
            i = i + 1
        i = index
        while(i < self._len):
            self._binaryGrayList[i] = tmp[i - index]
            i = i + 1
        '''
        i = index
        while(i < self._len):
            self._binaryGrayList[i], gene._binaryGrayList[i] = gene._binaryGrayList[i], self._binaryGrayList[i]  
            i = i + 1
    def Variation(self, index):
        self._binaryGrayList[index] = int(self._binaryGrayList[index]) ^ 1 
        
    def Equal(self, gene):
        for i in range(self._len):
            if(self._binaryGrayList[i] != gene._binaryGrayList[i]):
                return False
        return True

class Body:
    def __init__(self, limitList = [], precisionList = []):
        self._limitList = limitList
        self._precisionList = precisionList
        self._dimention = len(limitList)
        self._Set = [] # 记录支配的个体
        self._Np = 0 # 记录被支配数
        self._crow = 0
        self._level = -1

    def Init(self):
        self._DNA = []
        for i in range(self._dimention):
            g = Gene(self._limitList[i], self._precisionList[i])
            g.Random()
            self._DNA.append(g)    

    def Copy(self):
        b = Body(self._limitList, self._precisionList)
        b._DNA = []
        b._fitNess = self._fitNess
        for i in range(self._dimention):
            g = self._DNA[i].Copy()
            b._DNA.append(g)
        return b

    def Decode(self):
        values = []
        for i in range(self._dimention):
            values.append(self._DNA[i].Decode())
        return values
    
    def SetFitness(self, fitNess):
        self._fitNess = fitNess

    def GetFitNess(self):
        return self._fitNess

    def AddDomination(self, b):
        self._Set.append(b)

    def GetDomination(self):
        return self._Set

    def AddNp(self):
        self._Np = self._Np + 1

    def SubNp(self):
        self._Np = self._Np - 1

    def GetNp(self):
        return self._Np

    def SetLevel(self, level):
        self._level = level
    
    def GetLevel(self):
        return self._level

    def SetCrow(self, c):
        self._crow = c
    
    def GetCrow(self):
        return self._crow

    def ReSetSAndN(self):
        self._Np = 0
        self._Set = []
        self._level = -1

    def GetNegeCount(self):
        sum = 0
        for i in range(len(self._DNA)):
            sum += self._DNA[i].GetLen()
        return sum

    def CrossOver(self, body):
        #print("body cross Over")
        for i in range(self._dimention):
            if(random.random() > 0.9):
                continue
            len = self._DNA[i].GetLen()
            index = (int(random.random() * MAX)) % len
            self._DNA[i].CrossOver(body._DNA[i], index)

    def Variation(self, index):
        for i in range(len(self._DNA)):
            tmplen = self._DNA[i].GetLen()
            if (tmplen > index):
                self._DNA[i].Variation(index)
                break
            index -= tmplen

    def Equal(self, body):
        for i in range(self._dimention):
            if(not self._DNA[i].Equal(body._DNA[i])):
                return False
        return True

class NSGAII:
    def __init__(self, objectionFunctionList,limitList, precisionList, population, generation, needArchive = False):
        self._objectionFunctionList = objectionFunctionList
        self._limitList = limitList
        self._precisionList = precisionList
        self._population = population
        self._generation = generation
        self._lastGroup = [] # 上一代种群
        self._variationRate = 0.005
        self._paretoSet = []
        self._needArchive = needArchive
        self._archive = []

    def _getObjectionFunctionsValues(self, values):
        res = []
        for i in range(len(self._objectionFunctionList)):
            res.append((self._objectionFunctionList[i])(values))
        return res

    def _InitGroup(self): # 初始化种群
        self._group = []
        for _ in range(self._population):
            b = Body(self._limitList, self._precisionList)
            b.Init()
            self._group.append(b)

        self._lastGroup = []
        for _ in range(self._population):
            b = Body(self._limitList, self._precisionList)
            b.Init()
            values = b.Decode()
            values = self._getObjectionFunctionsValues(values)
            b.SetFitness(values)
            self._lastGroup.append(b)

    def _Domination_small(self, a, b):
        size = len(a)
        flag = False
        for i in range(size):
            if a[i] > b[i]:
                return False
            if a[i] < b[i]:
                flag = True
        return flag

    def _Domination_big(self, a, b):
        size = len(a)
        flag = False
        for i in range(size):
            if a[i] < b[i]:
                return False
            if a[i] > b[i]:
                flag = True
        return flag

    def _FastNondominatedSort(self, group, justFitst = False):
        # 将每个个体的S/N/Level清零
        for i in range(len(group)):
            group[i].ReSetSAndN()
        # 先计算每个个体支配的个体集合S 以及每个个体被支配的数量N
        size = len(group)
        H = []
        F = []
        level = 0
        for i in range(size):
            for j in range(size):
                if(i == j):
                    continue
                if(self._Domination_small(group[i].GetFitNess(), group[j].GetFitNess())):
                    group[i].AddDomination(group[j])
                elif(self._Domination_small(group[j].GetFitNess(), group[i].GetFitNess())):
                    group[i].AddNp()
                    if(justFitst):
                        break
            if(group[i].GetNp() == 0):
                group[i].SetLevel(level)
                H.append(group[i])
        if(justFitst):
            return [H]
        while(len(H) != 0):
            #print(len(H), ts)
            TMP = []
            for i in range(len(H)):
                TMP.append(H[i].Copy())
            F.append(TMP.copy())
            level = level + 1
            TMP = []
            for i in range(len(H)):
                tmpSet = H[i].GetDomination()
                for j in range(len(tmpSet)):
                    tmpSet[j].SubNp()
                    if(tmpSet[j].GetNp() == 0):
                        tmpSet[j].SetLevel(level)
                        TMP.append(tmpSet[j])
            H = TMP.copy()

        return F

    def _quickSortWithIndexValue(self, group, index, left, right):
        if(left >= right):
            return
        tmp = group[left].GetFitNess()[index]
        b = group[left]
        pl = left
        pr = right
        while(pl < pr):
            while(group[pr].GetFitNess()[index] > tmp and pr > pl):
                pr = pr - 1
            if(pr == pl):
                break
            group[pl] = group[pr]
            while(group[pl].GetFitNess()[index] <= tmp and pr > pl):
                pl = pl + 1
            if(pl == pr):
                break
            group[pr] = group[pl]
        group[pl] = b
        self._quickSortWithIndexValue(group, index, left, pl - 1)
        self._quickSortWithIndexValue(group, index, pr + 1, right)

    def _quickSortWithCrow(self, group, left, right):
        if(left >= right):
            return
        tmp = group[left].GetCrow()
        b = group[left]
        pl = left
        pr = right
        while(pl < pr):
            while(group[pr].GetCrow() > tmp and pr > pl):
                pr = pr - 1
            if(pr == pl):
                break
            group[pl] = group[pr]
            while(group[pl].GetCrow() <= tmp and pr > pl):
                pl = pl + 1
            if(pl == pr):
                break
            group[pr] = group[pl]
        group[pl] = b
        self._quickSortWithCrow(group, left, pl - 1)
        self._quickSortWithCrow(group, pr + 1, right)

    def _SortWithValueOnIndex(self, group, index):
        self._quickSortWithIndexValue(group, index, 0, len(group)-1)


    def _SortWithCrow(self, group):
        self._quickSortWithCrow(group, 0, len(group) - 1)


    def _CrowdingDistanceAssignment(self, group):
        # 计算每个个体的拥挤距离
        # 按照每个目标函数排序
        for i in range(len(group)):
            group[i].SetCrow(0)
        for i in range(len(self._objectionFunctionList)):
            self._SortWithValueOnIndex(group, i)
            for j in range(len(group)):
                if(j == 0):
                    group[j].SetCrow(999999999)
                elif(j == len(group) - 1):
                    t = group[j].GetCrow()
                    group[j].SetCrow(t + group[j].GetFitNess()[i] - group[j-1].GetFitNess()[i])
                else:
                    t = group[j].GetCrow()
                    group[j].SetCrow(t - group[j-1].GetFitNess()[i] + group[j+1].GetFitNess()[i])

        # 按照拥挤度排序
        self._SortWithCrow(group)
        #print(group[0]._DNA[0]._binaryGrayList)
        return group

    def _Copy(self): # 拷贝  先计算适应度，然后进行非支配排序，然后进行抽样拷贝
        #print("开始拷贝")

        for i in range(self._population):
            values = self._group[i].Decode()
            res = self._getObjectionFunctionsValues(values)
            self._group[i].SetFitness(res)
        tgroup = self._group + self._lastGroup
        # 把每个对象都拷贝一次 不然会出错
        group = []
        for i in range(len(tgroup)):
            group.append(tgroup[i].Copy())
        self._lastGroup = []
        for i in range(self._population):
            self._lastGroup.append(self._group[i].Copy())
        # 快速非支配排序
        F = self._FastNondominatedSort(group)
        # 拥挤度排序
        for i in range(len(F)):
            F[i] = self._CrowdingDistanceAssignment(F[i])
        count = 0
        tmpGroup = []

        ts = 0
        for i in range(len(F)):
            ts += len(F[i])
        if(ts != self._population * 2):
            print("出现错误{}_{}_{}".format(ts, self._population,len(group)))
            exit(0)

        for i in range(len(F)):
            for j in range(len(F[i])):
                tmpGroup.append(F[i][j].Copy())
                count = count + 1
                if(count == self._population):
                    break
            if(count == self._population):
                break
        self._group = tmpGroup
        
        for i in range(len(tmpGroup)):
            if(self._needArchive):
                if(not self._bodyInList(tmpGroup[i], self._archive)):
                    self._archive.append(tmpGroup[i].Copy())
            if(not self._bodyInList(tmpGroup[i], self._paretoSet)):
                self._paretoSet.append(tmpGroup[i].Copy())
        F = self._FastNondominatedSort(self._paretoSet, True)
        self._paretoSet = F[0]

        if(len(self._group) != self._population):
            print("出现错误")
            exit(1)

    def _bodyInList(self, body, List):
        for i in range(len(List)):
            if(body.Equal(List[i])):
                return True
        return False

    def _CrossOver(self): # 交叉
        #print("开始交叉")
        # 先两两配对
        half = int(self._population / 2)
        # 随机搭配
        selectedSet = set()
        for _ in range(half):
            a = (int(random.random() * MAX)) % self._population
            while(a in selectedSet):
                a = (a + 1) % self._population
            selectedSet.add(a)
            b = (int(random.random() * MAX)) % self._population
            while(b in selectedSet):
                b = (b + 1) % self._population
            selectedSet.add(b)
            self._group[a].CrossOver(self._group[b])       
    
    def _Variation(self): # 变异
        #print("开始变异")
        geneCount = self._group[0].GetNegeCount()
        totalCount = geneCount * self._population
        variationCount = int(totalCount * self._variationRate)
        for _ in range(variationCount):
            index = (int(random.random()) * MAX) % totalCount
            bodyNum = int(index / geneCount)
            index = index % geneCount
            #print("第{}个个体在{}发生变异".bodyNum, index)
            self._group[bodyNum].Variation(index)

    def Fit(self):
        self._InitGroup()
        for i in range(self._generation):
            print("开始第{}代".format(i))
            self._Copy()
            print("pareto解集中点的个数[{}]".format(len(self._paretoSet)))
            self._CrossOver()
            self._Variation()
        self._Copy()
        res = []
        for i in range(len(self._paretoSet)):
            res.append(self._paretoSet[i].Decode())
        if(self._needArchive):
            archive = []
            for i in range(len(self._archive)):
                archive.append(self._archive[i].Decode())
            return res, archive
        return res

def funcA(xList):
    Sum = 0
    for i in range(len(xList)):
        Sum += xList[i]**2
    return Sum

def funcB(xList):
    Sum = 0
    for i in range(len(xList)):
        Sum += (xList[i]-2)**2
    return Sum

if __name__ == "__main__":
    nsga = NSGAII([funcA, funcB], [[-10,10]], [10], 200, 30, True)
    g, a = nsga.Fit()
    print(g)
    X = []
    Y = []
    for i in range(len(g)):
        X.append(funcA(g[i]))
        Y.append(funcB(g[i]))

    plt.scatter(X,Y)
    plt.show()

    X = []
    Y = []
    for i in range(len(a)):
        X.append(funcA(a[i]))
        Y.append(funcB(a[i]))

    plt.scatter(X,Y)
    plt.show()