import numpy as np
import random

from numpy import linalg

from point import *
from face_new import *
from convert import *
from global_proporties import *
from vector_alg import *

class Boundary:
    def __init__(self,xmin,xmax,ymin,ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def flip_faces(self, faces):
        bx1 = False
        by1 = False
        bx2 = False
        by2 = False
        for f in faces:
            if bx1 and bx2 and by1 and by2:
                f.is_boundary_face = True