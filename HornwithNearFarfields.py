from __future__ import division
import meep as mp
print(mp.inf)
from matplotlib import pyplot as plt
import numpy as np
from meep.materials import Al
resolution = 10
sxy = 18
dpml = 1
cell = mp.Vector3(30,30,30)
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
thickness = 0.2
height = 1.4
zmin_bottom = -1.2
zmax_bottom = zmin_bottom + thickness
zmin_sides = zmin_bottom
zmax_sides = zmin_bottom + height
zmin_top = zmax_sides - thickness
zmax_top = zmax_sides
cell_zmin = zmin_bottom - height
cell_zmax = zmax_top + height
zmin_plug = zmin_bottom
zmax_plug = zmin_bottom + height
delta = mp.Vector3(-8,0,0)
bottom = mp.get_GDSII_prisms(mp.metal, gdsIIfile, bottom_layer, zmin_bottom, zmax_bottom)
sides = mp.get_GDSII_prisms(mp.metal, gdsIIfile, side_layer, zmin_sides, zmax_sides)
top = mp.get_GDSII_prisms(mp.metal, gdsIIfile, top_layer, zmin_top, zmax_top)
plug = mp.get_GDSII_prisms(mp.metal, gdsIIfile, plug_layer, zmin_plug, zmax_plug)
bottom[0].center += delta
top[0].center += delta
sides[0].center += delta
sides[1].center += delta
plug[0].center += delta
geometry = sides + plug + bottom + top
sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=df),center=mp.Vector3(-7,0,-0.6),component=src_cmpt,size=mp.Vector3(0.0,0.2,0.0))]
sim = mp.Simulation(resolution=resolution,cell_size=cell,boundary_layers=pml_layers,geometry=geometry,sources=sources)
nearfield_box = sim.add_near2far(
    fcen,
    df,
    64,
    mp.Near2FarRegion(center=mp.Vector3(0,+0.5*sxy*0.66,0),size=mp.Vector3(sxy,0,0), direction=mp.X),
    mp.Near2FarRegion(center=mp.Vector3(0,-0.5*sxy*0.66,0),size=mp.Vector3(sxy,0,0), direction=mp.X, weight=-1),
    mp.Near2FarRegion(center=mp.Vector3(+0.5*sxy,0,0),size=mp.Vector3(0,sxy*0.66,0), direction=mp.Y),
    mp.Near2FarRegion(center=mp.Vector3(-0.5*sxy,0,0), size=mp.Vector3(0,sxy*0.66,0), direction=mp.Y,weight=-1),
    mp.Near2FarRegion(center=mp.Vector3(0,0,+0.5*sxy), size=mp.Vector3(0,0,sxy), direction=mp.Z),
    mp.Near2FarRegion(center=mp.Vector3(0,0,-0.5*sxy), size=mp.Vector3(0,0,sxy), direction=mp.Z, weight=-1))
sim.run(mp.to_appended("ey",mp.at_every(0.5,mp.output_efield_y)),until=200)
# fig1 = plt.figure(dpi=150)
# sim.plot2D(ax=fig1.gca(),eps_parameters={'cmap': 'gray'})
# plt.show()
which_freq = 10
which_real_freq = 10.0
npts = 360
npts2 = 180
angles = 2 * np.pi/npts * np.arange(npts)
# print(len(angles))
theta = np.pi/npts * np.arange(npts)
# print(theta)
# print(len(theta))
r = 1000
E = np.zeros((npts, 3), dtype=np.complex128)
H = np.zeros((npts, 3), dtype=np.complex128)
freq_marker_1 = (which_freq - 1) * 6
freq_marker_2 = freq_marker_1 + 3
# for n in range(npts):
#     ff = sim.get_farfield(nearfield_box, mp.Vector3(r * np.cos(angles[n]), r * np.sin(angles[n])))
#     E[n, :] = [np.conj(ff[j]) for j in range(freq_marker_1, freq_marker_2)]
#     H[n, :] = [ff[j + 3] for j in range(freq_marker_1, freq_marker_2)]

for n in range(npts):
    ff = sim.get_farfield(nearfield_box, mp.Vector3(r*np.sin(theta[n]), 0, r * np.cos(theta[n])))
    E[n, :] = [np.conj(ff[j]) for j in range(freq_marker_1, freq_marker_2)]
    H[n, :] = [ff[j + 3] for j in range(freq_marker_1, freq_marker_2)]

Px = np.real(E[:, 1] * H[:, 2] - E[:, 2] * H[:, 1])
Py = np.real(E[:, 2] * H[:, 0] - E[:, 0] * H[:, 2])
Pz = np.real(E[:, 0] * H[:, 1] - E[:, 1] * H[:, 0])
Pr = np.sqrt(np.square(Px) + np.square(Py) + np.square(Pz))
directivity = 10.0 * np.log10(Pr/max(Pr)) #decibels
print(E)
print(H)
# print(Px)
# print(Py)
print(Pz)
# print(Pr)
# print(directivity)
fig3 = plt.figure("Radiation Pattern")
# plt.plot(angles, directivity, color='blue')
plt.plot(theta, directivity, color='red')
# ax = plt.gca()
# ax.set_rlim(-30, 0)
# ax.set_rticks([-30, -15, -3, 0])
# ax.grid(True)
# ax.set_rlabel_position(180)
# ax.tick_params(labelsize=18)
# ax.set_title("Frequency (GHz): " + "{:3.1f}".format(which_real_freq))
# ax.title.set_fontsize(18)
plt.xlabel('theta(radians)')
plt.ylabel('dB')
plt.title('Radiation Pattern at Frequency (GHz): ' + '{:3.1f}'.format(which_real_freq))
plt.show()
