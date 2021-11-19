import numpy as np

from point import *

class Face1D:
    def __init__(self, number, center:Point) -> None:
        self.number = number
        self.center = center

    def __str__(self):
        return f'Face number {self.number} between %s'%(self.center)