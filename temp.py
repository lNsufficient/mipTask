import numpy
import scipy.linalg as sl
import matplotlib.cm as cm
import matplotlib.pyplot as plt


class Room(object):
    def __init__(self, divisor):
        self.defaultValues()
        self.dx = 1/divisor
        self.columns = divisor 
        self.rows = self.columns - 1 
        self.values = numpy.ones(self.columns*self.rows)*self.t0
        self.bvStart = self.boundaryValuesSetup()*divisor**2
        self.omega = self.createOmega(numpy.size(self.values))*divisor**2

    def defaultValues(self):    
        self.t0 = 15
        self.heat = 40
        self.window = 5
        self.ohm = 0.8
    
    def ohmUpdate(self):
        self.values = self.values*self.ohm + (1-self.ohm)*self.oldValues


    def derivative(self,u1,u2):
        return (u1-u2)/self.dx

    def nextSolution(self, bvNew):
        #print("omega : ", self.omega)
        self.values, self.oldValues = sl.solve(self.omega, self.newBV(bvNew)), self.values
    
    def newBV(self, bvNew):
        bv = self.bvStart + self.computeBoundaryValues(bvNew)
##        bv = numpy.ones(numpy.size(self.values))*(-self.t0)
        #print("BV used:", bv)
        return bv

#TODO: punkterna med oklara boundaryValues har inte fått sina värden. 

class Room2(Room):
    def __init__(self, divisor):
        self.defaultValues()
        self.dx = 1/divisor
        self.columns = divisor - 1 #antal element som beräknas i varje rad
        self.rows = self.columns*2 + 1
        self.values = numpy.ones(self.columns*self.rows)*self.t0
        self.bvStart =self.boundaryValuesSetup()/(self.dx**2)
        self.omega = self.createOmega(numpy.size(self.values))/(self.dx**2)
       # print("bvStart :", self.bvStart)

    def newBV(self, bv1to2, bv3to2):
    #    print("bv1to2: ", bv1to2)
     #   print("bv3to2: ", bv3to2)
        bv = self.bvStart + self.computeBoundaryValues(bv1to2, bv3to2)/(self.dx**2)
  #      bv = numpy.ones(numpy.size(self.values))*(-self.t0)
      #  print("Bv used (room2):", bv)
        return bv
   
    def nextSolution(self, bv1to2, bv3to2):
        self.values, self.oldValues = sl.solve(self.omega, self.newBV(bv1to2, bv3to2)), self.values
        #print("Temperatures in room 2: ", self.values) 

    def sendBoundaryValues(self):
        bv2to1 = numpy.zeros(self.columns)
        bv2to3 = numpy.zeros(self.columns)
       
        for i in range(self.columns):
            index1 = (self.columns+1)*self.columns+i*self.columns
            index3 = self.columns - 1 + i*self.columns
            bv2to1[i] = self.values[index1]
            bv2to3[i] = self.values[index3]
  
        #print("bv till 1: ", bv2to1)
        #print("bv till 3: ", bv2to3)
        return bv2to1, bv2to3

    def createOmega(self,size):
        diag = -4*numpy.ones(size)
        subdiag1 = numpy.ones(size-1)
        for i in range(size-1):
            if (i+1)%self.columns == 0:
                subdiag1[i]= 0
        subdiag2 = numpy.ones(size-self.columns)
        Omega2 = numpy.diag(diag) + numpy.diag(subdiag1,1) + numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,self.columns) +  numpy.diag(subdiag2,-self.columns)
        return Omega2

    def boundaryValuesSetup(self):
            
        bv = numpy.zeros(numpy.size(self.values))
        bv[:self.columns] -= self.heat
        bv[-self.columns:] -= self.window
        for i in range(self.columns+1): 
            bv[i*self.columns] -= self.t0
            bv[self.columns*(self.columns+1) + i*self.columns -1] -= self.t0
        return bv

    def computeBoundaryValues(self, bv1to2, bv3to2):
        bv = numpy.zeros(numpy.size(self.bvStart))
        for i in range(self.columns):
            bv[self.columns*(i+1)-1] -= bv3to2[i]

        for i in range(self.columns):
            bv[self.columns*(self.columns+1)+self.columns*i] -= bv1to2[i]
        
        return bv
        

class Room1(Room):
   
    def boundaryValuesSetup(self):
        bv = numpy.zeros(numpy.size(self.values))
        bv[:self.columns] -= self.t0
        bv[-self.columns:] -= self.t0
        for i in range(self.rows):
            bv[i*self.columns] -= self.heat
        return bv
            
    def createOmega(self, size):
        diag = -4*numpy.ones(size)
        subdiag1 = numpy.ones(size-1)
        for i in range(size-1):
            if (i+1)%self.columns ==0:
                subdiag1[i] = 0
        subdiag2 = numpy.ones(size-self.columns)
        for i in range(size):
            if i%self.columns == self.columns-1:
                diag[i] +=1
        Omega1 = numpy.diag(diag) + numpy.diag(subdiag1,1)+ numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,self.columns)+ numpy.diag(subdiag2,-self.columns)
        return Omega1
        

    def computeBoundaryValues(self, bv2to1):
        bv = numpy.zeros(numpy.size(self.bvStart))
        for i in range(self.rows):
            #bv2to1 contains values, derivative is calculated in this room.
            index = self.columns*(i+1) - 1
            bv[index] = -1/self.dx*self.derivative(bv2to1[i],self.values[index])
            #self.bv[self.columns*(i+1) - 1] = -1/(2*dx)*(self.values[self.columns-1 + self.columns*i] - bv2to1[i])
        return bv

    def sendBoundaryValues(self):
        bv1to2 = numpy.zeros(self.rows)
        for i in range(self.rows):
            index = (i+1)*self.columns-1
            bv1to2[i] = self.values[index]
        return bv1to2
            
class Room3(Room):

    def boundaryValuesSetup(self):
        bv = numpy.zeros(numpy.size(self.values))
        bv[:self.columns] -= self.t0
        bv[-self.columns:] -= self.t0
        for i in range(self.rows):
            bv[i*self.columns + self.columns - 1] -= self.heat
        return bv

    def computeBoundaryValues(self, bv2to3):
        bv = numpy.zeros(numpy.size(self.bvStart))
        for i in range(self.rows):
            index = self.columns*i
            bv[index] = -self.derivative(bv2to3[i], self.values[index])/self.dx
        return bv   
    
    def createOmega(self, size):
        diag = -4*numpy.ones(size)
        subdiag1 = numpy.ones(size-1)
        for i in range(size-1):
            if (i+1)%self.columns == 0:
                subdiag1[i] = 0
        for i in range(size):
            if i%self.columns ==0:
                diag[i] += 1
        subdiag2 = numpy.ones(size-self.columns)
        Omega3 = numpy.diag(diag) + numpy.diag(subdiag1,1)+ numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,self.columns)+ numpy.diag(subdiag2,-self.columns)
        print(Omega3)
        return Omega3   
    
    def sendBoundaryValues(self):
        bv3to2 = numpy.zeros(self.rows)
        for i in range(self.rows):
            index = self.columns * i
            bv3to2[i] = self.values[index]
        return bv3to2

class Apartment(object):
    def __init__(self, divisor):
        self.r1 = Room1(divisor)
        self.r2 = Room2(divisor)
        self.r3 = Room3(divisor)

    def __call__(self, iterations):
        bv1to2, bv3to2 = self.r2.t0*numpy.ones(self.r2.columns), self.r2.t0*numpy.ones(self.r2.columns) 
        for i in range(iterations):
            self.r2.nextSolution(bv1to2, bv3to2)
            bv2to1, bv2to3 = self.r2.sendBoundaryValues()
            self.r1.nextSolution(bv2to1)
            self.r3.nextSolution(bv2to3)
            self.r1.ohmUpdate()
            self.r2.ohmUpdate()
            self.r3.ohmUpdate()
            bv1to2, bv3to2 = self.r1.sendBoundaryValues(), self.r3.sendBoundaryValues()
        print(self.fullApartment())
        return self.r1.values, self.r2.values, self.r3.values


    def fullApartment(self):
        values = numpy.ones((self.r2.rows + 2, self.r2.columns*3 + 4))*self.r2.t0
        values[0,self.r2.columns+1:(self.r2.columns)*2 + 3] = self.r2.heat
        values[-1,self.r2.columns+2:(self.r2.columns)*2 + 3] = self.r2.window
        values[:self.r3.rows+2,-1] = self.r3.heat
        values[-(self.r1.rows+2):,0] = self.r1.heat
        for i in range(self.r1.rows):
            values[self.r1.rows+2+i,1:1+self.r1.columns] = self.r1.values[i*self.r1.columns:self.r1.columns*(i+1)]
        for i in range(self.r2.rows):
            values[1+i, self.r2.columns+2:self.r2.columns*2+2] = self.r2.values[i*self.r2.columns:self.r2.columns*(i+1)]
        for i in range(self.r3.rows):
            values[1+i, self.r2.columns*2+2:-1] = self.r3.values[i*self.r3.columns:self.r3.columns*(i+1)]
        plt.imshow(values, extent=(0, 3, 0, 2), interpolation='nearest')
        plt.colorbar()
        plt.show()
       
        return values

        

    
apt = Apartment(20)
r1, r2, r3 = apt(20)   

print("r1: ", r1)
print("r2: ", r2)
print("r3: ", r3)
