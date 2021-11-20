from convert import *
from global_proporties import *
from point import *
from cell1d import *
from face1d import *
from convert1d import *

def create_mesh(n=30):
    dim = 1
    points = np.linspace(0,dim,n)

    points_obj = []
    faces = []
    cells = []

    for i, p in enumerate(points):
        points_obj.append(Point(np.array([p,0])))
        faces.append(Face1D(i, Point(np.array([p,0]))))
    
    for i in range(len(points_obj)-1):
        new_cell = Cell(i)
        new_cell.boundary_points.append(points_obj[i])
        new_cell.boundary_points.append(points_obj[i+1])
        cells.append(new_cell)

    for i, c in enumerate(cells):
        c.neighbors.append(c[i-1])
        c.neighbors.append(c[i+1])

    return cells, faces

if __name__=='__main__':
    cells, faces = create_mesh(n=30)
    
    for c in cells:
        print(c)