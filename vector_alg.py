import numpy as np
import matplotlib.pyplot as plt
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def display_vector(center, vector, color='k', scale=40):
    plt.quiver(*center, vector[0], vector[1], width=0.001, 
                headlength=1, headwidth=3, scale=scale, color=color)  

def get_tangent(x,y,dxdy, width=0.1):
    x_tan = np.linspace(x-width,x+width,10)
    print(x_tan)
    y_tan = dxdy*(x_tan - np.ones([1,10])[0]*x) + np.ones([1, 10])[0]*y
    print(y_tan)
    return [x_tan, y_tan]

def plot_tangent(cell):
    

    tan = get_tangent(cell.center.X, cell.rho, cell.gradients[0], width=1/8 * cell.longest_side)
   
    plt.plot(tan[0], tan[1], 'r')
