from numpy import true_divide
from point import *

class Face:
    def __init__(self, p1:Point, p2:Point):
        if p1 is p2:
            raise Exception("Face generation with 2 equal points")

        # geometry
        self.center = Point()
        self.boundary_points = [p1, p2]
        self.calc_center()
        self.surface = p1.distance(p2)

        # primitive values
        self.rho = 0
        self.p = 0
        self.u = 0
        self.v = 0
    
    
    def calc_center(self):
        boundary_points_ko =np.zeros([len(self.boundary_points),2])
        k = 0
        for p in self.boundary_points:
            boundary_points_ko[k] = p.getValue()
            k+=1
        self.center  = boundary_points_ko.mean(axis=0)

    def is_equal_to(self, other_face):
        # Takes another Face() as input and checks its points
        # if the points are equal the face declares itself as equal
        i = 0 # i = 0 not identical; i = 1 one point identical; i = 2 both points identical 
        p1 = self.boundary_points[0]
        p2 = self.boundary_points[1]
        p1o = other_face.boundary_points[0]
        p2o = other_face.boundary_points[1]
        
        if p1 == p1o:
            i +=1
        if p1 == p2o:
            i +=1
        if p2 == p1o:
            i +=1
        if p2 == p2o:
            i +=1
        
        if i == 2:
            return True
        else:
            return False
        

    def __str__(self):
        return f"Face between {self.boundary_points[0].__str__()} and {self.boundary_points[1].__str__()} center : {self.center}"

if __name__ == "__main__":
    p1 = Point(np.array([1,1]))
    p2 = Point(np.array([3,2]))
    f1 = Face(p1,p2)
    f2 = Face(p2,p1)
    print(f1.is_equal_to(f2))