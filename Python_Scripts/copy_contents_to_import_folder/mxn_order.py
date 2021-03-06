import math

class Arithmetic:

    def length(self,N):
        return N

    def lengthstr(self,N):
        return "%.1f"%self.length(N)

    def clusters(self,N):
        tol=10**-9
        maxy=int(N)
        res=[]
        for y in range(2,maxy+1):
            x=int(2*N)-y #int needed for decimal orders
            xf=float(2*N)-y
            diff=xf-x
            if diff < tol:
                if x >= y:
                    res.append((x,y))
        
        return res

class Geometric:

    def length(self,N):
        return math.sqrt(N)

    def lengthstr(self,N):
        return "%.4f"%self.length(N)

    def clusters(self,N):
        tol=10**-9
        maxy=int(N/2)
        res=[]
        for y in range(2,maxy+1):
            x=N/y
            xf=float(N)/y
            diff=xf-x
            if diff < tol:
                if x >= y:
                    res.append((x,y))

        return res
        
class Quadratic:

    def length(self,N):
        return math.sqrt(0.5*N)

    def lengthstr(self,N):
        return "%.4f"%self.length(N)

    def clusters(self,N):
        from math import sqrt
        tol=10**-9
        maxy=int(sqrt(N))
        res=[]
        for y in range(2,maxy+1):
            x2=N-y*y
            x=sqrt(x2)
            xint=int(x)
            diff=x-xint
            if diff < tol:
                if x >= y:
                    res.append((xint,y))

        return res
