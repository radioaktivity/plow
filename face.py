from point import *

class Face:
    def __init__(self, p1:Point, p2:Point):
        # geometry
        self.center = 0
        self.boundary_points = [p1, p2]
        self.sruface = 0

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

    def __str__(self):
        return f"Face between {self.boundary_points[0].__str__()} and {self.boundary_points[1].__str__()} center : {self.center}"

if __name__ == "__main__":
    p1 = Point(np.array([1,1]))
    p2 = Point(np.array([3,2]))
    f1 = Face(p1,p2)
    f1.calc_center()
    print(f1)