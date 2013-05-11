from geometry import *
from solver import *

geometry = Geometry()
solver = Solver(geometry)

solver.convergeScalarFlux(order=8, num_mesh=100, tol=1E-3, max_iters=50)
