from geometry import *
from solver import *
import numpy as np

geometry = Geometry()
solver = Solver(geometry)

meshes = [250]
orders = [24]

dancoffs = np.zeros((len(orders), len(meshes)))
flux_ratios = np.zeros((len(orders), len(meshes)))

for i in range(len(orders)):
    for j in range(len(meshes)):
        print 'order = %d, mesh = %d' % (orders[i], meshes[j])
        solver.computeDancoff(meshes[j], orders[i])
        dancoffs[i,j] = solver.getDancoffFactor()
        flux_ratios[i,j] = solver.getFluxRatio()
        solver.plotAngularFlux()
