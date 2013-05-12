import math
import numpy
import matplotlib.pyplot as plt


class Geometry(object):

    def __init__(self, pitch=1.26, radius=0.4):
        self.pitch = pitch
        self.radius = radius
        self.fuel_width = math.sqrt(math.pi * self.radius**2)
        self.createFuelBox()


    ############################################################################
    ###########################  GETTERS / SETTERS  ############################
    ############################################################################

    def setPitch(self, pitch):
        self.pitch = pitch
        self.createFuelBox()

    def setRadius(self, radius):
        self.radius = radius
        # Compute an equivalent rectangular fuel width

        self.fuel_width = math.sqrt(math.pi * self.radius**2)

        self.createFuelBox()

    def getPitch(self):
        return self.pitch

    def getRadius(self):
        return self.radius

    def getFuelWidth(self):
        return self.fuel_width

    def getFuelArea(self):
        return self.fuel_width**2

    def getModeratorArea(self):
        return self.pitch**2 - self.fuel_width**2

    def getNumMesh(self):
        return self.num_mesh

    def getMeshCellDelta(self):
        return self.delta

    def getMeshStencil(self):
       return self.mesh_stencil

    def getMeshCellArea(self):
        return self.cell_area

    def getXIndex(self, x):
        return math.floor(x / self.delta) - 1

    def getYIndex(self, y):
        return math.floor(y / self.delta) - 1

    ############################################################################
    ############################  MESH GENERATION  #############################
    ############################################################################

    def createFuelBox(self):
            
        pitch_minus_radius = self.pitch - self.fuel_width

        self.fuel_min_x = pitch_minus_radius * 0.5
        self.fuel_min_y = pitch_minus_radius * 0.5
        self.fuel_max_x = pitch_minus_radius * 0.5 + self.fuel_width
        self.fuel_max_y = pitch_minus_radius * 0.5 + self.fuel_width
        
        self.fuel_top_left = numpy.array([self.fuel_min_x, self.fuel_max_y])
        self.fuel_top_right = numpy.array([self.fuel_max_x, self.fuel_max_y])
        self.fuel_bottom_left = numpy.array([self.fuel_min_x, self.fuel_min_y])
        self.fuel_bottom_right = numpy.array([self.fuel_max_x, self.fuel_min_y])


    def fuelContainsCoords(self, coords):
        
        if not isinstance(coords, list):
            coords = [coords]

        contains = []

        for i in range(len(coords)):

            contains.append(True)

            if coords[i][0] > self.fuel_max_x:
                contains[i] = False
            if coords[i][0] < self.fuel_min_x:
                contains[i] = False
            if coords[i][1] > self.fuel_max_y:
                contains[i] = False
            if coords[i][1] < self.fuel_min_y:
                contains[i] = False
                
        return contains


    def fuelContainsMeshCell(self, mesh_cell):
        return all(self.fuelContainsCoords(mesh_cell))


    def fuelOverlapsMeshCell(self, mesh_cell):
        return any(self.fuelContainsCoords(mesh_cell))


    def computeFuelMeshCellArea(self, mesh_cell):

        if self.fuelContainsMeshCell(mesh_cell):
            return self.cell_area
        else:
            width_x = 0.
            width_y = 0.

            # Fuel contains top left corner of mesh cell
            if self.fuelContainsCoords(mesh_cell[0])[0]:

                # Fuel contains top left and top right corners of mesh cell
                if self.fuelContainsCoords(mesh_cell[1])[0]:
                    width_x = self.delta
                    width_y = mesh_cell[0][1] - self.fuel_min_y
                # Fuel contains top left and bottom left corners of mesh cell
                elif self.fuelContainsCoords(mesh_cell[3])[0]:
                    width_y = self.delta
                    width_x = self.fuel_max_x - mesh_cell[0][0]
                # Fuel contains top left corner only of mesh cell
                else:
                    width_x = self.fuel_max_x - mesh_cell[0][0]
                    width_y = mesh_cell[0][1] - self.fuel_min_y

            # Fuel contains top right corner of mesh cell
            elif self.fuelContainsCoords(mesh_cell[1])[0]:

                # Fuel contains top right and top left corners of mesh cell
                if self.fuelContainsCoords(mesh_cell[0])[0]:
                    width_x = self.delta
                    width_y = mesh_cell[1][1] - self.fuel_min_y
                # Fuel contains top right and bottom right corners of mesh cell
                elif self.fuelContainsCoords(mesh_cell[2])[0]:
                    width_y = self.delta
                    width_x = mesh_cell[1][0] - self.fuel_min_x
                # Fuel contains top right corner only of mesh cell
                else:
                    width_x = mesh_cell[1][0] - self.fuel_min_x
                    width_y = mesh_cell[1][1] - self.fuel_min_y

            # Fuel contains bottom right corner of mesh cell
            elif self.fuelContainsCoords(mesh_cell[2])[0]:
                # Fuel contains bottom right and top right corners of mesh cell
                if self.fuelContainsCoords(mesh_cell[1])[0]:
                    width_x = mesh_cell[2][0] - self.fuel_min_x
                    width_y = self.delta
                # Fuel contains bottom right and bottom left corners of mesh cell
                elif self.fuelContainsCoords(mesh_cell[3])[0]:
                    width_x = self.delta
                    width_y = self.fuel_max_y - mesh_cell[2][1]
                # Fuel contains bottom right corner only of mesh cell
                else:
                    width_x = mesh_cell[2][0] - self.fuel_min_x
                    width_y = self.fuel_max_y - mesh_cell[2][1]

            # Fuel contains bottom left corner of mesh cell
            elif self.fuelContainsCoords(mesh_cell[3])[0]:
                # Fuel contains bottom left and top left corners of mesh cell
                if self.fuelContainsCoords(mesh_cell[0])[0]:
                    width_x = self.fuel_max_x - mesh_cell[3][0]
                    width_y = self.delta
                # Fuel contains bottom left and bottom right corners of mesh cell
                elif self.fuelContainsCoords(mesh_cell[2])[0]:
                    width_x = self.delta
                    width_y = self.fuel_max_y - mesh_cell[3][1]
                # Fuel contains bottom left corner only of mesh cell
                else:
                    width_x = self.fuel_max_x - mesh_cell[3][0]
                    width_y = self.fuel_max_y - mesh_cell[3][1]

            return width_x * width_y

        

    def generateMeshStencil(self, num_mesh=15):

        self.num_mesh = num_mesh
        self.mesh_stencil = numpy.zeros((self.num_mesh, self.num_mesh))
        self.delta = self.pitch / self.num_mesh
        self.cell_area = self.delta**2

        fuel_area = 0.
                
        # Create array of #'s in [0,1] with ratio of moderator (0) to fuel (1)
        for i in range(num_mesh):
            for j in range(num_mesh):
                
                # Create this mesh cell
                top_left = numpy.array([self.delta * i, self.delta * (j+1)])
                top_right = numpy.array([self.delta * (i+1), self.delta * (j+1)])
                bottom_right = numpy.array([self.delta * (i+1), self.delta * j])
                bottom_left = numpy.array([self.delta * i, self.delta * j])
                mesh_cell = [top_left, top_right, bottom_right, bottom_left]
                
                # Compute area bounded by fuel
                if not self.fuelOverlapsMeshCell(mesh_cell):
                    fuel_area = 0.
                else:
                    fuel_area = self.computeFuelMeshCellArea(mesh_cell)
                    
                # compute area ratio
                if fuel_area != 0.:
                    self.mesh_stencil[i][j] = fuel_area / self.cell_area


    #############################################################################
    ################################  PLOTTING  #################################
    #############################################################################

    def plotMeshStencil(self):

        plt.figure()

        plt.pcolor(numpy.linspace(0, self.pitch, self.num_mesh+1),
                   numpy.linspace(0, self.pitch, self.num_mesh+1),
                   self.mesh_stencil, edgecolors='k', linewidth=0.25)

        plt.axis([0, self.pitch, 0, self.pitch])
        plt.title('Ordinates Pin Cell Mesh')
        plt.savefig('pin-cell-mesh-' + str(self.num_mesh) + '.png')
