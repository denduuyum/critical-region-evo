
class UnionFind():
    def __init__(self, n):
        self.n = n
        self.parents = [-1] * n
        self.sz = [1] * n

    def find(self, x):
        if self.parents[x] < 0:
            return x
        else:
            self.parents[x] = self.find(self.parents[x])
            return self.parents[x]

    def union(self, x, y):
        x = self.find(x)
        y = self.find(y)

        if x == y:
            return

        if self.sz[x] > self.sz[y]:
            x, y = y, x

        self.sz[x] += self.sz[y]
        self.parents[y] = x

    def size(self, x):
        return self.sz[self.find(x)]

    def same(self, x, y):
        return self.find(x) == self.find(y)

    def sizes(self):
        isCalc = [0] * self.n
        csize = []
        for i in range(self.n):
            x = self.find(i)
            # print('uf: ', i, x)

            if isCalc[x] == 0:
                isCalc[x] = 1
                csize.append(self.size(x))
        return csize

    def eval(self):
        t = 0
        isCalc = [0] * self.n
        for i in range(self.n):
            x = self.find(i)
            if isCalc[x] == 0:
                isCalc[x] = 1
                
                t += self.size(x) * (self.size(x) - 1) / 2

        return t
                
