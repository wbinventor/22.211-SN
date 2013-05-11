import geometry
import quadrature
import numpy as np
import matplotlib.pyplot as plt
import math


class Solver(object):

    def __init__(self, geometry):

        self.geometry = geometry
        self.quadrature = quadrature.LevelSymmetricQuadrature()
        self.fuel_total_xs = 100.
        self.moderator_total_xs = 0.25


    def generateQuadratureSet(self):

        self.quad_set = self.quadrature.getQuadratureSet(self.quad_order)
        self.na = self.quad_set['num angles per octant']
        self.wgt = self.quad_set['weight']


    def initializeMeshCellSources(self):

        mesh_cell_area = self.geometry.getMeshCellArea()
        fuel_area = self.geometry.getFuelArea()
        self.source = self.mesh_stencil * (mesh_cell_area / fuel_area)
        self.source /= (4. * math.pi)


    def initializeMeshCellXSData(self):

        # Initialize arrays for the homogenized mesh cell cross-section data
        self.total_xs = self.mesh_stencil * self.fuel_total_xs
        self.total_xs += (1.0 - self.mesh_stencil) * self.moderator_total_xs


    def initializeAngularFluxes(self):

        # Indices: cell edge(0 - left, 1 - bottom, etc), mesh index, angle
        self.wave = np.zeros((4, self.num_mesh, self.na*2))

        # Indices: quadrant, x, y, angle
        self.ang_flux = np.zeros((self.num_mesh, self.num_mesh, 4, self.na))


    def initializeScalarFluxes(self):

        self.scalar_flux = np.zeros((self.num_mesh, self.num_mesh))
        self.old_scalar_flux = np.ones((self.num_mesh, self.num_mesh))

    def computeScalarFlux(self):

        self.scalar_flux.fill(0.)

        self.scalar_flux = self.wgt[np.newaxis, np.newaxis, np.newaxis,...] * \
                           self.ang_flux

#        for y in self.mesh:
#            for x in self.mesh:

                
                # Loop over angles for this mesh cell
#                for quad in np.arange(4):
#                    for ang in np.arange(self.na):
                    
#                        self.scalar_flux[y,x] += self.wgt[ang] * \
#                                                self.ang_flux[y,x,quad,ang]


    def convergeScalarFlux(self, num_mesh=15, order=4, max_iters=250, tol=1E-3):
        
        # Initialize mesh cell xs data and sources
        self.num_mesh = num_mesh
        self.mesh = np.arange(self.num_mesh)

        self.geometry.generateMeshStencil(self.num_mesh)
        self.delta = self.geometry.getMeshCellDelta()
        self.mesh_stencil = self.geometry.getMeshStencil()
        self.initializeMeshCellSources()
        self.initializeMeshCellXSData()

        # Get a level-symmetric quadrature set for this order        
        self.quad_order = order
        self.generateQuadratureSet()

        # Pre-compute mu and eta fractions
        self.mus = 2. * self.quad_set['mu'] / self.delta
        self.etas = 2. * self.quad_set['eta'] / self.delta

        # Initialize angular and scalar flux arrays
        self.initializeAngularFluxes()
        self.initializeScalarFluxes()

        # Initialize the initial residual on the cell-avg scalar fluxes
        eps = np.inf

        # Loop over transport sweep until converged
        for i in range(max_iters):

            print 'Iteration %d with eps = %f' % (i, eps)
            
            # Plot the scalar flux from this iteration
            self.plotScalarFlux(iteration=i)
            
            # Perform a transport sweep
            self.transportSweep()
            
            # Compute scalar flux from angular fluxes
            self.computeScalarFlux()

            # Normalize the scalar flux to the average
            self.scalar_flux /= np.average(self.scalar_flux.flatten())

            # Check for convergence
            eps = np.sqrt(np.sum((self.scalar_flux - self.old_scalar_flux)**2))
                          
            # If converged, break the loop
            if (eps < tol):
                break
            # Otherwise, store old scalar flux and repeat
            else:
                self.old_scalar_flux = np.copy(self.scalar_flux)

        print 'converged in %d iterations with eps = %f' % (i, eps)



    def transportSweep(self, reflect=True):

        # Zero the angular flux
        self.ang_flux.fill(0.)

        quadrants = np.arange(4)

        for quad in quadrants:
            self.sweepQuadrant(quad, reflect)


    def sweepQuadrant(self, quad, reflect=True):

        if quad == 0:
            step_x = +1
            step_y = +1
            side_x = 1
            side_y = 0
            offset_x = 0
            offset_y = 0
        elif quad == 1:
            step_x = -1
            step_y = -1
            side_x = 3
            side_y = 2
            offset_x = self.na
            offset_y = self.na
        elif quad == 2:
            step_x = -1
            step_y = +1
            side_x = 1
            side_y = 2
            offset_x = self.na
            offset_y = 0
        else:
            step_x = +1
            step_y = -1
            side_x = 3
            side_y = 0
            offset_x = 0
            offset_y = self.na

        # Loop over mesh cells
        for y in self.mesh[::step_y]:
            for x in self.mesh[::step_x]:

                # Get the total xs and source for this mesh cell
                total_xs = self.total_xs[y][x]
                source = self.source[y][x]

                # Loop over angles for this mesh cell
                for ang in np.arange(self.na):

                    # Get mu and eta for this angle
                    mu = self.mus[ang]
                    eta = self.etas[ang]

                    # offset is for wavefront
                    # no need for offset for ang_flux since quad does it
                    ang_offset_x = ang + offset_x
                    ang_offset_y = ang + offset_y

                    # Compute cell-averaged angular flux (diamond difference)
                    self.ang_flux[y,x,quad,ang] = source
                    self.ang_flux[y,x,quad,ang] += \
                                             mu*self.wave[side_y,y,ang_offset_y]
                    self.ang_flux[y,x,quad,ang] += \
                                             eta*self.wave[side_x,x,ang_offset_x]
                    self.ang_flux[y,x,quad,ang] /= (total_xs + mu + eta)

                    # Compute flux on right/left cell edge of next cell
                    self.wave[side_x,x,ang_offset_x] = \
                                             2.*self.ang_flux[y,x,quad,ang] - \
                                             self.wave[side_x,x,ang_offset_x]

                    # Compute flux on top/bottom cell edge of next cell
                    self.wave[side_y,y,ang_offset_y] = \
                                             2.*self.ang_flux[y,x,quad,ang] - \
                                             self.wave[side_y,y,ang_offset_y]

        # Reflective boundary conditions
        if reflect:
            for i in self.mesh:
                for ang in np.arange(self.na):
                    self.wave[side_x-2,i,ang_offset_x] = \
                                             self.wave[side_x,i,ang_offset_x]
                    self.wave[side_y-2,i,ang_offset_y] = \
                                             self.wave[side_y,i,ang_offset_y]


    def plotScalarFlux(self, iteration):

        plt.figure()
        plt.pcolor(np.linspace(0, self.geometry.pitch, self.num_mesh+1),
                   np.linspace(0, self.geometry.pitch, self.num_mesh+1),
                   self.scalar_flux)        
        plt.axis([0, self.geometry.getPitch(), 0, self.geometry.getPitch()])
        
        plt.title('SN Pin Cell Flux (iteration ' + str(iteration)  +')')

        suffix = '-order-' + str(self.quad_order) + '-iter-' + str(iteration)
        plt.savefig('flux' + suffix + '.png')
