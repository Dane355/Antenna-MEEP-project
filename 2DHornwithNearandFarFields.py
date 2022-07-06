from __future__ import division
import meep as mp
from matplotlib import pyplot as plt
resolution = 10
sxy = 18
dpml = 1
cell = mp.Vector3(30,30)
pml_layers = [mp.PML(dpml)]
fcen = 0.1
df = 0.075
src_cmpt = mp.Ey
gdsIIfile = 'horn_layout2.gds'
bottom_layer = 1
side_layer = 2
top_layer = 3
cell_layer = 4
plug_layer = 5
epsilon_horn = mp.metal
sides = mp.get_GDSII_prisms(epsilon_horn, gdsIIfile, side_layer)
plug = mp.get_GDSII_prisms(epsilon_horn, gdsIIfile, plug_layer)
delta = mp.Vector3(-8,0)
sides[0].center += delta
sides[1].center += delta
plug[0].center += delta
geometry = sides + plug
sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=df), center=mp.Vector3(-7,0,0), component=src_cmpt, size=mp.Vector3(0.0,0.2,0.0))]
sim = mp.Simulation(resolution=resolution, cell_size=cell, boundary_layers=pml_layers, geometry=geometry, sources=sources)
nearfield_box = sim.add_near2far(fcen, df, 64, mp.Near2FarRegion(center=mp.Vector3(0,+0.5*sxy*0.66), size=mp.Vector3(sxy,0)), mp.Near2FarRegion(center=mp.Vector3(0,-0.5*sxy*0.66), size=mp.Vector3(sxy,0), weight=-1), mp.Near2FarRegion(center=mp.Vector3(+0.5*sxy,0), size=mp.Vector3(0,sxy*0.66)), mp.Near2FarRegion(center=mp.Vector3(-0.5*sxy,0), size=mp.Vector3(0,sxy*0.66), weight=-1))
sim.run(mp.to_appended("ey",mp.at_every(0.5,mp.output_efield_y)),until=200)
fig1 = plt.figure(dpi=150)
sim.plot2D(ax=fig1.gca(),eps_parameters={'cmap':'gray'})
plt.show()