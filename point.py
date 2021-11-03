import math
import numpy as np



class Point(object):
    '''Creates a point on a coordinate plane with values x and y.'''

    COUNT = 0

    def __init__(self, *arg):
        '''Defines x and y variables'''    
        if arg == ():
            self.X = None
            self.Y = None
        elif type(arg[0]) is np.ndarray:
            self.X = arg[0][0]
            self.Y = arg[0][1]
        else: 
            raise Exception("Wrong Argument for Point() was transmitted")

    def move(self, dx, dy):
        '''Determines where x and y move'''
        self.X = self.X + dx
        self.Y = self.Y + dy

    def __str__(self):
        return "Point(%s,%s)"%(self.X, self.Y) 

    def getValue(self):
        return(np.array([self.X, self.Y]))

    def getX(self):
        return self.X

    def getY(self):
        return self.Y

    def distance(self, other):
        dx = self.X - other.X
        dy = self.Y - other.Y
        return math.sqrt(dx**2 + dy**2)
    

if __name__ == "__main__":
    p1 = Point(np.array([0,0]))
    p2 = Point(np.array([0,0]))
    z = p1
    del z
    print(p1)
    print(p2)
    print(p1 == z)
