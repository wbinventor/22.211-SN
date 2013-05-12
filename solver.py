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


    ############################################################################
    ###########################  GETTERS / SETTERS #############################
    ############################################################################

    def getDancoffFactor(self):
        return self.dancoff


    def getIsolatedFuelReactionRate(self):
        return self.isolated_RR

    
    def getLatticeFuelReactionRate(self):
        return self.lattice_RR


    def getFluxRatio(self):
        return self.flux_ratio


    ############################################################################
    ############################  INITIALIZATION  ##############################
    ############################################################################

    def initializeSolver(self, num_mesh=15, order=4, max_iters=250, tol=1E-3):

        print 'initialize solver'

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


    ############################################################################
    ############################  SOLVER ROUTINES  #############################
    ############################################################################

    def computeDancoff(self, num_mesh=15, order=4, max_iters=250, tol=1E-3):

        self.initializeSolver(num_mesh, order, max_iters, tol)

        # Perform a transport sweep using vacuum boundary conditions to compute
        # the isolated reaction rate in the fuel
        self.transportSweep(reflect=False)
        self.computeScalarFlux()
#        self.plotScalarFlux(suffix='vacuum')
        self.isolated_RR = self.computeFuelReactionRate()

        # Compute scalar flux from angular fluxes
        self.convergeScalarFlux(num_mesh, order, max_iters, tol)
        self.computeScalarFlux()
        self.total_RR = self.computeTotalReactionRate()
        self.lattice_RR = self.computeFuelReactionRate()
        self.flux_ratio = self.computeFluxRatio()

        self.dancoff = 1. - (1. - self.lattice_RR / self.total_RR) / \
                            (1. - self.isolated_RR / self.total_RR)

        print 'isolated RR = %1.10f, total RR = %1.10f, lattice RR = %1.10f' % \
               (self.isolated_RR, self.total_RR, self.lattice_RR)
        print 'dancoff = %f, flux ratio = %f' % (self.dancoff, self.flux_ratio)
        
        self.dancoff = 1. - (4.*math.pi*

    def convergeScalarFlux(self, num_mesh=15, order=4, max_iters=250, tol=1E-3):

        self.initializeSolver(num_mesh, order, max_iters, tol)
        
        # Initialize the initial residual on the cell-avg scalar fluxes
        eps = np.inf

        # Loop over transport sweep until converged
        for i in range(max_iters):

            print 'Iteration %d with eps = %f' % (i, eps)
                        
            # Perform a transport sweep
            self.transportSweep()
            
            # Compute scalar flux from angular fluxes
            self.computeScalarFlux()
            self.normalizeScalarFlux()

            # Plot the scalar flux from this iteration
            self.plotScalarFlux(iteration=i)

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


    ############################################################################
    ############################  DATA PROCESSING  #############################
    ############################################################################

    def computeScalarFlux(self):

        self.scalar_flux = (self.ang_flux * self.wgt[np.newaxis, np.newaxis, 
                                                     np.newaxis, :])
        self.scalar_flux = self.scalar_flux.sum(axis=3).sum(axis=2)

    
    def normalizeScalarFlux(self):

        # Normalize the scalar flux to the average
        self.scalar_flux /= np.average(self.scalar_flux.flatten())


    def computeTotalReactionRate(self):

        total_RR = 0.

        for y in self.mesh:
            for x in self.mesh:
                total_RR += self.scalar_flux[y,x] * self.total_xs[y,x]

        total_RR *= self.geometry.getMeshCellArea()
        return total_RR


    def computeFuelReactionRate(self):

        fuel_RR = 0.

        for y in self.mesh:
            for x in self.mesh:
                if self.mesh_stencil[y,x] > 0.:
                    fuel_RR += self.scalar_flux[y,x] * self.total_xs[y,x]

        fuel_RR *= self.geometry.getMeshCellArea()

        return fuel_RR


    def computeFluxRatio(self):
        
        fuel_flux = np.average(self.scalar_flux, weights=self.mesh_stencil)
        moderator_flux = np.average(self.scalar_flux, 
                                    weights=(1.-self.mesh_stencil))

        return fuel_flux / moderator_flux



    ############################################################################
    ###############################  PLOTTING  #################################
    ############################################################################

    def plotScalarFlux(self, iteration=0, suffix=''):

        plt.figure()
        plt.pcolor(np.linspace(0, self.geometry.pitch, self.num_mesh+1),
                   np.linspace(0, self.geometry.pitch, self.num_mesh+1),
                   self.scalar_flux)        
        plt.axis([0, self.geometry.getPitch(), 0, self.geometry.getPitch()])
        
        plt.title('SN Pin Cell Flux (order = %d)' % (self.quad_order))

        suffix = '-order-%d-iter-%d' %(self.quad_order, iteration) + suffix
        plt.savefig('flux' + suffix + '.png')


    def plotAngularFlux(self):
        
        # Get indices for mu/eta at the polar angle (xi) nearest the xy plane
        min_xi = min(self.quad_set['xi'])
        min_xi_indices = []

        for index, value in enumerate(self.quad_set['xi']):
            if abs(value - min_xi) < 1E-5:
                min_xi_indices.append(index)

        # Compute coordinates for the locations of interest (1-5)
        pitch = self.geometry.getPitch()
        fuel_width = self.geometry.getFuelWidth()

        print 'pitch = %f, fuel width = %f' % (pitch, fuel_width)

        coords = [] # Stored as (y,x) tuples

        coords.append((pitch/2., pitch/2.))
        coords.append((pitch/2. + fuel_width/2., pitch/2. + fuel_width/2.))
        coords.append((pitch, pitch))
        coords.append((pitch/2., pitch/2. + fuel_width/2.)) 
        coords.append((pitch/2., pitch))

        # Compute the index into the mesh for the locations of interest
        coords_indices = []

        for coord in coords:
            indices = (self.geometry.getYIndex(coord[0]),
                       self.geometry.getXIndex(coord[1]))
            coords_indices.append(indices)

            
        # Extract and plot the angular flux for each coordinate
        plt.figure()

        for indices in coords_indices:

            y = indices[0]
            x = indices[1] 

            ang_fluxes = []
            thetas = []
            
            # Second quadrant (quadrant 1)
            for index in min_xi_indices:
                ang_fluxes.append(self.ang_flux[y,x,1,index])
                thetas.append(math.pi - math.acos(self.quad_set['mu'][index]))

            # First quadrant (quadrant 0)
            for index in min_xi_indices:
                ang_fluxes.append(self.ang_flux[y,x,0,index])
                thetas.append(math.acos(self.quad_set['mu'][index]))

            # Third quadrant(quadrant 2)
            for index in min_xi_indices:
                ang_fluxes.append(self.ang_flux[y,x,2,index])
                thetas.append(math.pi + math.acos(self.quad_set['mu'][index]))
            
            # Fourth quadrant(quadrant 3)
            for index in min_xi_indices:
                ang_fluxes.append(self.ang_flux[y,x,3,index])
                thetas.append(2*math.pi - math.acos(self.quad_set['mu'][index]))

            # Convert arrays to numpy arrays
            ang_fluxes = np.array(ang_fluxes)
            thetas = np.array(thetas)

            # Sort the thetas and fluxes in order of increasing theta
            theta_indices = np.array(np.argsort(thetas))

            thetas = thetas[theta_indices]
            ang_fluxes = ang_fluxes[theta_indices]

            # Normalize the fluxes to the average
            ang_fluxes /= np.average(ang_fluxes)

            plt.plot(thetas, ang_fluxes, linewidth=2)

        plt.xlabel('Theta [Radians]')
        plt.ylabel('Normalized Angular Flux')
        plt.legend(['Location 1', 'Location 2', 'Location 3', 
                    'Location 4', 'Location 5'])
        plt.title('S' + str(self.quad_order) + ' Pin Cell Angular Flux')
        plt.grid()
        plt.ylim([0.5, 2.])
        plt.xlim([0, 2.*math.pi])
        plt.show()
        suffix = '-order-%d-' % (self.quad_order)
        plt.savefig('angular-flux' + suffix + '.png')


