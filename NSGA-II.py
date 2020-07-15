# coding=utf-8
import random
## 单个个体
# 对于单个个体而言，每个个体有多个维度，每个维度有一定的取值区间，对于精度也有一定的要求
# 另外 每个个体应当保存自己的适应度

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

class Body:
    def __init__(self, limitList = [], precisionList = []):
        self._limitList = limitList
        self._precisionList = precisionList
        self._dimention = len(limitList)
        self._Set = [] # 记录支配的个体
        self._Np = 0 # 记录被支配数
        self._crow = 0

    def Init(self):
        self._DNA = []
        for i in range(self._dimention):
            g = Gene(self._limitList[i], self._precisionList[i])
            g.Random()
            self._DNA.append(g)    

    def Copy(self):
        b = Body(self._limitList, self._precisionList)
        b._DNA = []
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

class NSGAII:
    def __init__(self, objectionFunctionList,limitList, precisionList, population, generation):
        self._objectionFunctionList = objectionFunctionList
        self._limitList = limitList
        self._precisionList = precisionList
        self._population = population
        self._generation = generation
        self._lastGroup = [] # 上一代种群

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
            self._lastGroup.append(b)

    def _Domination(self, a, b):
        size = len(a)
        flag = False
        for i in range(size):
            if a[i] > b[i]:
                return False
            if a[i] < b[i]:
                flag = True
        return flag

    def _FastNondominatedSort(self, group):
        # 将每个个体的S和N清零
        for i in range(len(group)):
            group[i].ReSetSAndN()
        # 先计算每个个体支配的个体集合S 以及每个个体被支配的数量N
        size = len(group)
        for i in range(size):
            for j in range(size):
                if(i == j):
                    continue
                if(self._Domination(group[i].GetFitNess(), group[j].GetFitNess())):
                    group[i].AddDomination(group[j])
                elif(self._Domination(group[j].GetFitNess(), group[i].GetFitNess())):
                    group[i].AddNp()

        H = []
        F = []
        level = 0
        for i in range(size):
            if(group[i].GetNp() == 0):
                H.append(group[i])
                group[i].SetLevel(level)

        while(True):
            if(len(H) == 0):
                break
            F.append(H)
            level = level + 1
            TMP = []
            for i in range(len(H)):
                tmpSet = H[i].GetDomination()
                for j in range(len(tmpSet)):
                    tmpSet[j].SubNp()
                    if(tmpSet[j].GetNp() == 0):
                        tmpSet[j].SetLevel(level)
                        TMP.append(tmpSet[j])
            H = TMP
        '''
        for i in range(len(F)):
            for j in range(len(F[i])):
                print(F[i][j].GetFitNess())
            print("===",i)
        '''
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
        print(group[0]._DNA[0]._binaryGrayList)
        return group

    def _Copy(self): # 拷贝  先计算适应度，然后进行非支配排序，然后进行抽样拷贝
        for i in range(self._population):
            values = self._group[i].Decode()
            res = self._getObjectionFunctionsValues(values)
            self._group[i].SetFitness(res)
        group = self._group + self._lastGroup
        self._lastGroup = self._group
        # 快速非支配排序
        F = self._FastNondominatedSort(group)
        # 拥挤度排序
        for i in range(len(F)):
            F[i] = self._CrowdingDistanceAssignment(F[i])
        count = 0
        tmpGroup = []
        for i in range(len(F)):
            for j in range(len(F[i])):
                tmpGroup.append(F[i][j])
                count = count + 1
                if(count == self._population):
                    break
            if(count == self._population):
                break
        self._group = tmpGroup

    def _CrossOver(self): # 交叉
        a = 1

    def _Variation(self): # 变异
        a = 1

    def Fit(self):
        self._InitGroup()
        self._Copy()

def funcA(xList):
    Sum = 0
    for i in range(len(xList)):
        Sum += xList[i]
    return Sum

def funcB(xList):
    Sum = 0
    for i in range(len(xList)):
        Sum += xList[i] / (i + 1)
    return Sum

if __name__ == "__main__":
    nsga = NSGAII([funcA, funcB], [[0,10],[0,5]], [6,2], 10, 10)
    nsga.Fit()