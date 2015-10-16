class Room(Object):
    self.t0 = 15
    self.heat = 40
    self.window = 5
        
    def __init__(self, divisor):
        self.dx = 1/divisor
        self.columns = divisor 
        self.rows = columns - 1 
        self.bvStart = self.boundaryValuesSetup()
        self.values = zeroes 
        self.omega = self.createOmega(numpy.size(self.bv))

    def derivative(self,u1,u2,u3):
        return -1/
#TODO: punkterna med oklara boundaryValues har inte fått sina värden. 

class Room2(Room):
    def __init__(self, divisor):
        dx = 1/divisor
        columns = divisor - 1 #antal element som beräknas i varje rad
        rows = collumns*2 + 1
  
	def createOmega(size):
        diag = -4*numpy.ones(size)
        subdiag1 = numpy.array([ 1-i%2 for i in range(size-1)])
        subdiag2 = numpy.ones(size-2)
        Omega2 = numpy.diag(diag) + numpy.diag(subdiag1,1) + numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,2) +  numpy.diag(subdiag2,-2)
        return Omega2

    def boundaryValuesSetup()
            
        bv = numpy.zeroes(numpy.size(self.values))
        bv[:self.columns] -= self.heat
        bv[-self.columns:] -= self.window
        for i in range(columns+1): 
            bv[i*self.columns] -= self.t0
            bv[self.columns*(self.columns+1) + i*self.columns -1] -= self.t0
        return bv

    def updateBoundaryValues():
        self.bv = numpy.zeros(numpy.size(self.bvStart))
        for i in range(self.columns):
            self.bv[self.columns*(i+1)-1] -= bv3to2[i]

        for i in range(self.columns):
            self.bv[self.columns*(self.columns+1)+self.columns*i] -= bv1to2[i]
        

class Room1(Room):
   
    def boundaryValuesSetup():
        bv = numpy.zeroes(numpy.size(self.values))
        bv[:self.colums] -= self.t0
        bv[-self.columns:] -= self.t0
        for i in range(self.rows):
            bv[i*self.columns] -= self.heat
        return bv
            
    def createOmega(size):
        diag = -4*numpy.ones(size)
        subdiag1 = numpy.ones(size-1)
        subdiag2 = numpy.ones(size-3)
        diag[size/2-1] +=1
        diag[-1] +=1
        subdiag1[(size-2)/2] = 0
        Omega1 = numpy.diag(diag) + numpy.diag(subdiag1,1)+ numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,3)+ numpy.diag(subdiag2,-3)
        return Omega1
        

    def updateBoundaryValues():
        self.bv = numpy.zeros(numpy.size(self.bvStart))
        for i in range(rows):
            self.bv[self.columns*(i+1) - 1] = -1/(2*dx)*(self.values[self.columns-1 + self.columns*i] - bv2to1[i])
    
class Room3(Room):

    def boundaryValuesSetup():
        bv = numpy.zeroes(numpy.size(self.values))
        bv[:self.columns] -= self.t0
        bv[-self.columns:] -= self.t0
        for i in range(self.rows):
            bv[i*self.columns + self.columns - 1] -= self.heat
        return bv

    def updateBoundaryValues():
        self.bv = numpy.zeros(numpy.size(self.bvStart))
        for i in range(rows):
            self.bv[self.columns*i+1] = -self.derivative(bv2to3[i], self.values[self.columns*i])
       
    def createOmega(size):
        diag = -4*numpy.ones(size)
        subdiag1 = numpy.ones(size-1)
        subdiag2 = numpy.ones(size-3)
        diag[size/2] +=1
        diag[0] +=1
        subdiag1[(size-2)/2] = 0
        Omega3 = numpy.diag(diag) + numpy.diag(subdiag1,1)+ numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,3)+ numpy.diag(subdiag2,-3)
        return Omega3   
    
