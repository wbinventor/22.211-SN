from geometry import *
import math

test = Geometry()

print 'pitch = ' + str(test.getPitch())
print 'radius = ' + str(test.getRadius())
print 'fuel_width = ' + str(test.getFuelWidth())

test.generateMeshStencil(200)
print test.getMeshStencil()
test.plotMeshStencil()
