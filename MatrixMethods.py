import numpy

#Task 1

deltax = 1/3

def createOmega2(size):
    diag = -4*numpy.ones(size)
    subdiag1 = numpy.array([ 1-i%2 for i in range(size-1)])
    subdiag2 = numpy.ones(size-2)
    Omega2 = numpy.diag(diag) + numpy.diag(subdiag1,1) + numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,2) +  numpy.diag(subdiag2,-2)
    return Omega2

size = 10
size2 = 6

#Omega1 = 1/(deltax**2)*numpy.array([[-4,1,1.,0],[1,-3,0,1],[1,0,-4,1],[0,1,1,-3]])
#Omega2 = createOmega2(size)
#Omega3 = 1/(deltax**2)*numpy.array([[-3,1,1.,0],[1,-4,0,1],[1,0,-3,1],[0,1,1,-4]])


def createOmega1(size):
    diag = -4*numpy.ones(size)
    subdiag1 = numpy.ones(size-1)
    subdiag2 = numpy.ones(size-3)
    diag[size/2-1] +=1
    diag[-1] +=1
    subdiag1[(size-2)/2] = 0
    Omega1 = numpy.diag(diag) + numpy.diag(subdiag1,1)+ numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,3)+ numpy.diag(subdiag2,-3)
    return Omega1
    
def createOmega3(size):
    diag = -4*numpy.ones(size)
    subdiag1 = numpy.ones(size-1)
    subdiag2 = numpy.ones(size-3)
    diag[size/2] +=1
    diag[0] +=1
    subdiag1[(size-2)/2] = 0
    Omega3 = numpy.diag(diag) + numpy.diag(subdiag1,1)+ numpy.diag(subdiag1,-1) + numpy.diag(subdiag2,3)+ numpy.diag(subdiag2,-3)
    return Omega3

    
Omega1 = createOmega1(size2)
Omega2 = createOmega2(size)
Omega3 = createOmega3(size2)    

print(Omega1)
print(Omega2)
print(Omega3)

######################

#Constant right hand sides

b1 = numpy.array([-55,-15,-15,-55,-15,-10])
b2 = numpy.array([-27.5,-15,-55,-15,-15,-55])
b3 = numpy.array([-55,-40,-15,0,-15,-15,0,-15,-5,-20])

#More general
def createb1(row):
    b1 = numpy.zeros(row*(row-1))
    b1[:row] +=15.
    b1[-row:-1] +=15
    for i in range(row*(row-1)):
        if (i+1)%row==1:
            b1[i] +=40
    b1[-1] = 10
    return b1
    
def createb3(row):
    b3 = numpy.zeros(row*(row-1))
    b3[0]=27.5
    b3[1:row] +=15
    b3[-row:] +=15
    for i in range(row*(row-1)):
        if (i+1)%row==0:
            b3[i] +=40
    return b3
b1 = createb1(3)
b3 = createb3(7)
print(b1)
print(b3)
