'''''''''''''''''''''
# @FileName:EMOMMO.py
# @author:ZhaoXinYi
# @version:0.0.1
# @Date:2020.07.17
# @BSD
'''''''''''''''''''''
import NSGAII_FOREMOMMO as nsga
import PeekDetection
import JADE
import math

class EMOMMO:
    def __init__(self, objectionFunction, limit):
        self._obejectionFunction = objectionFunction
        self._limit = limit
        size = len(limit)
        self._precisionList = []
        for _ in range(size):
            self._precisionList.append(5)
    
    # 求最大化问题
    def Fit(self):
        ga = nsga.NSGAII_FOREMOMMO(self._obejectionFunction, self._limit, self._precisionList, 100, 50, True)
        _, fitNessLandscape = ga.Fit()
        vlist = []
        for i in range(len(fitNessLandscape)):
            vlist.append(self._obejectionFunction(fitNessLandscape[i]))
        
        peekSet = PeekDetection.PeekDection(fitNessLandscape, vlist)
        print("从[{}]个点中，找到[{}]个峰值集合".format(len(fitNessLandscape), len(peekSet)))
        bestPoint = []
        for i in range(len(peekSet)):
            max = -9999999999
            bestPoint.append(peekSet[i][0])
            for j in range(len(peekSet[i])):
                if(self._obejectionFunction(peekSet[i][j]) >= max):
                    bestPoint[i] = peekSet[i][j]
                    max = self._obejectionFunction(peekSet[i][j])
        localRangeList = []
        for i in range(len(bestPoint)):
            localRange = []
            for j in range(len((self._limit))):
                space = self._limit[j][1] -  self._limit[j][0]
                left = bestPoint[i][j] - space * 0.025
                right = bestPoint[i][j] + space * 0.025
                if(left < self._limit[j][0]):
                    left = self._limit[j][0]
                if(right > self._limit[j][1]):
                    right = self._limit[j][1]
                localRange.append([left, right])
            localRangeList.append(localRange)
        solution = []
        values = []

        def maxToMin(x, function = self._obejectionFunction):
            return -function(x)

        for i in range(len(localRangeList)):
            print("在[{}]内进行localSearch".format(localRangeList[i]))
            j = JADE.MyJADE(maxToMin, localRangeList[i], 30, 50)
            s, b = j.Fit()
            solution.append(s)
            values.append(-b)
        return solution, values

def testObjection(x):
    return math.sin(x[0]*4) + math.cos(x[1]*4)

if __name__ == "__main__":
    e = EMOMMO(testObjection, [[0,10],[0,10]])
    s, b = e.Fit()
    for i in range(len(s)):
        print("第{}个最优解[{}] 函数值[{}]".format(i,s[i], b[i]))
    