# by Andrew Householder
import meep as mp
import numpy as np

# can define cell in this script or in gdsII file
# cell_size = mp.Vector3(10,10,10)
pml_layers = [mp.PML(0.25)]

gdsIIfile = 'horn_layout.gds'
bottom_layer = 1
side_layer = 2
top_layer = 3
cell_layer = 4

thickness = 0.2
height = 0.95  # same as a

zmin_bottom = -1.2
zmax_bottom = zmin_bottom + thickness
zmin_sides = zmin_bottom
zmax_sides = zmin_bottom + height
zmin_top = zmax_sides - thickness
zmax_top = zmax_sides

cell_zmin = zmin_bottom - height
cell_zmax = zmax_top + height

# dielectric constant of horn
epsilon_horn = mp.metal

# defining the meep geoemtric objects and their properties
cell = mp.GDSII_vol(gdsIIfile, cell_layer, cell_zmin, cell_zmax)
bottom = mp.get_GDSII_prisms(epsilon_horn, gdsIIfile, bottom_layer, zmin_bottom, zmax_bottom)
sides = mp.get_GDSII_prisms(epsilon_horn, gdsIIfile, side_layer, zmin_sides, zmax_sides)
top = mp.get_GDSII_prisms(epsilon_horn, gdsIIfile, top_layer, zmin_top, zmax_top)

# geometry is a list of meep geoemtric objects
geometry = bottom + sides + top

# visulaization script -> source not included
sources = [mp.Source(mp.CustomSource(src_func=0),
                     component=mp.Ex,
                     center=mp.Vector3(0, 0, 0),
                     size=mp.Vector3(0, 1, 0))]

sim = mp.Simulation(resolution=20,
                    cell_size=cell.size,
                    boundary_layers=pml_layers,
                    geometry=geometry)

sim.init_sim()

eps_data = sim.get_epsilon()

# getting the epsilon (dielectric constant) data of the objects and
# saving in numpy array .npy file
np.save('horn_data.npy', eps_data)

'''
from mayavi import mlab
s = mlab.contour3d(eps_data, colormap="YlGnBu")
mlab.show()
'''