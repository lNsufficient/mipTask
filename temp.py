class room(Object):
    def __init__(self, divisor):
        dx = 1/divisor
        columns = divisor 
        rows = columns - 1

    def derivative(self,u1,u2,u3):
        return -1/


class Room2:
    def __init__(self, divisor):
        dx = 1/divisor
        columns = divisor - 1 #antal element som beräknas i varje rad
        rows = collumns*2 + 1
  

    def boundaryValuesRoom2():
            
        bv = numpy.zeroes(numpy.size(self.values))
        bv[:self.columns] = -self.heat
        bv[-self.columns:] = -self.window
        for i in range(columns+1): 
            bv[i*self.columns] = -self.t0
            bv[self.columns*(self.columns+1) + i*self.columns -1] = -self.t0
        self.bv = bv

    def updateBoundaryValuesRoom2():
        for i in range(self.columns):
            self.bv[self.columns*(i+1)-1] = -bv3to2[i]

        for i in range(self.columns):
            self.bv[self.columns*(self.columns+1)+self.columns*i] = -bv1to2[i]
        

class room1:
   
    def boundaryValuesRoom1():
        bv = numpy.zeroes(numpy.size(self.values))
        bv[:self.colums] = -self.t0
        bv[-self.columns:] = -self.t0
        for i in range(self.rows):
            bv[i*self.columns] = -self.heat
        

    def updateBoundaryValuesRoom1():
        for i in range(rows):
            self.bv[self.columns*(i+1) - 1] = -1/(2*dx)*(self.values[self.columns-1 + self.columns*i] - bv2to1[i])
    
class room3:

    def boundaryValuesRoom3():
        bv = numpy.zeroes(numpy.size(self.values))
        bv[:self.columns] = -self.t0
        bv[-self.columns:] = -self.t0
        for i in range(self.rows):
            bv[i*self.columns + self.columns - 1] = -self.heat
        

    def updateBoundaryValuesRoom3():
        for i in range(rows):
            self.bv[self.columns*i+1] = -self.derivative(bv2to3[i], self.values[self.columns*i])
       
    
