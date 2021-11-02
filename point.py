import math
import numpy as np
class Point(object):
    '''Creates a point on a coordinate plane with values x and y.'''

    COUNT = 0

    def __init__(self, array):
        '''Defines x and y variables'''
        self.X = array[0]
        self.Y = array[1]

        # not relevant value
        self.U = np.zeros([2,1])
        self.p = 0

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