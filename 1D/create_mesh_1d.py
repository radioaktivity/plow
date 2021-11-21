from convert import *
from global_proporties import *
from point import *
from cell1d import *
from face1d import *
from convert1d import *

def create_mesh(n=30):
    dim = 100
    points = np.linspace(0,dim,n)

    points_obj = []
    faces = []
    cells = []

    for i, p in enumerate(points):
        points_obj.append(Point(np.array([p,0])))
        faces.append(Face1D(i, Point(np.array([p,0]))))
    
    faces[0].wormhole_face = faces[-1]
    faces[-1].wormhole_face = faces[0]

    for i in range(len(points_obj)-1):
        new_cell = Cell(i)
        new_cell.boundary_points.append(points_obj[i])
        new_cell.boundary_points.append(points_obj[i+1])
        cells.append(new_cell)
    
    for i,c in enumerate(cells):
        c.faces.append(faces[i])
        c.faces.append(faces[i+1])

    for i, c in zip(range(len(cells)-1), cells):
        c.neighbors.append(cells[i-1])
        c.neighbors.append(cells[i+1])
    cells[-1].neighbors.append(cells[-2])
    cells[-1].neighbors.append(cells[0])

    for c in cells:
        c.calc_all()


    return cells, faces

if __name__=='__main__':
    cells, faces = create_mesh(n=30)
    

    for c in cells:
        print(c)
    for f in faces:
        print(f)