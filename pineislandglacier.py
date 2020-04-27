# %%
import os
import sys

sys.path.append(f"{os.environ['ISSM_DIR']}/bin")
sys.path.append(f"{os.environ['ISSM_DIR']}/lib")
sys.path.append(f"{os.environ['ISSM_DIR']}/src/m/dev")

# assert os.environ["ISSM_DIR"] == "/opt/issm/trunk"
os.environ["PYTHONSTARTUP"] = f"{os.environ['ISSM_DIR']}/src/m/dev/devpath.py"

import devpath

# import ISSM

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pygmt
import xarray as xr

# # %matplotlib inline

from model import model
from bamg import bamg
from export_netCDF import export_netCDF
from loadmodel import loadmodel
from plotmodel import plotmodel
from InterpFromGridToMesh import InterpFromGridToMesh
from savevars import savevars
from parameterize import parameterize
from setflowequation import setflowequation
from verbose import verbose
from toolkits import toolkits
from generic import generic
from solve import solve

from cuffey import cuffey

if not os.path.exists("Models"):
    os.makedirs(name="Models")

# %%
# Step 1 - Mesh Generation
domain = "pig_domain.exp"
assert os.path.exists(domain)
hinit = 10000  # element size for the initial mesh (originally 10000)
hmax = 40000  # maximum element size of the final mesh (originally 40000)
hmin = 5000  # minimum element size of the final mesh (originally 5000)
gradation = 1.7  # maximum size ratio between two neighboring elements
err = 8  # maximum error between interpolated and control field

# Generate an initial uniform mesh (resolution = hinit m)
md = bamg(model(), "domain", domain, "hmax", hinit)

# Load velocities
with xr.open_dataset(
    filename_or_obj="Data/antarctic_ice_vel_phase_map_v01.nc"
) as nsidc_vel:
    x = nsidc_vel.x.data
    y = np.flipud(m=nsidc_vel.y.data)
    vx = np.flipud(m=nsidc_vel.VX.data.astype(np.float64))
    vy = np.flipud(m=nsidc_vel.VY.data.astype(np.float64))

# Interpolate velocities onto coarse mesh
vx_obs = InterpFromGridToMesh(
    x=x, y=y, data=vx, x_mesh=md.mesh.x, y_mesh=md.mesh.y, default_value=0
)
vy_obs = InterpFromGridToMesh(
    x=x, y=y, data=vy, x_mesh=md.mesh.x, y_mesh=md.mesh.y, default_value=0
)
vel_obs = np.sqrt(vx_obs[0] ** 2 + vy_obs[0] ** 2)

# Adapt the mesh to minimize error in velocity interpolation
# https://issm.jpl.nasa.gov/documentation/tutorials/mesh
md = bamg(
    md,
    "hmax",
    hmax,
    "hmin",
    hmin,
    "gradation",
    gradation,
    "field",
    vel_obs,
    "err",
    err,
)

# ploting
plotmodel(md, "data", "mesh")

# Save model
export_netCDF(md=md, filename="Models/pig_mesh_generation.nc")

# %%
# Step 2 - Masks
md = loadmodel(path="Models/pig_mesh_generation.nc")

# Load ALBMAP dataset
# https://doi.pangaea.de/10.1594/PANGAEA.734145
with xr.open_dataset(filename_or_obj="Data/ALBMAPv1.nc") as albmap:
    x1 = albmap.x1.data.astype(np.float64)
    y1 = albmap.y1.data.astype(np.float64)
    mask = albmap.mask.data.astype(np.float64)

# interpolate onto our mesh vertices
groundedice = InterpFromGridToMesh(x1, y1, mask, md.mesh.x, md.mesh.y, 0)[0]
groundedice[groundedice <= 0] = -1

# fill in the md.mask structure
md.mask.groundedice_levelset = groundedice  # ice is grounded for mask equal one
md.mask.ice_levelset = -1 * np.ones(
    shape=md.mesh.numberofvertices
)  # ice present when negative

# plotting
plotmodel(
    md,
    "data",
    md.mask.groundedice_levelset,
    "title",
    "grounded/floating",
    # "data",
    # md.mask.ice_levelset,
    # "title",
    # "ice/no-ice",
)

# Save model
export_netCDF(md=md, filename="Models/pig_setmask.nc")

# %%
# Step 3 - Parameterization
md = loadmodel(path="Models/pig_setmask.nc")
md = parameterize(md=md, parametername="Pig/Pig_py.par")

# Save model
export_netCDF(md=md, filename="Models/pig_parameterization.nc")

# %%
# Step 4 - Control Method inverting for basal friction
md = loadmodel(path="Models/pig_parameterization.nc")
md.miscellaneous.name = "pig_control_drag_fullstokes"  # give model a name!

# Extrude Mesh
md = model.extrude(md, 9, 3)  # 9 layers, with extrusion exponent of 3
assert md.mesh.numberoflayers == 9

# Control general
md.inversion.iscontrol = 1
md.inversion.maxsteps = 20
md.inversion.maxiter = 40
md.inversion.dxmin = 0.1
md.inversion.gttol = 1.0e-4
md.verbose = verbose("control", True)

# Cost functions
md.inversion.cost_functions = [101, 103, 501]
md.inversion.cost_functions_coefficients = np.ones(shape=(md.mesh.numberofvertices, 3))
md.inversion.cost_functions_coefficients[:, 0] = 1
md.inversion.cost_functions_coefficients[:, 1] = 1
md.inversion.cost_functions_coefficients[:, 2] = 8e-15

# Controls
md.inversion.control_parameters = ["FrictionCoefficient"]
md.inversion.min_parameters = 1 * np.ones(shape=(md.mesh.numberofvertices, 1))
md.inversion.max_parameters = 200 * np.ones(shape=(md.mesh.numberofvertices, 1))

# Additional parameters
md.stressbalance.restol = 0.01
md.stressbalance.reltol = 0.1
md.stressbalance.abstol = np.NaN

# Set Flow Equation
# https://issm.jpl.nasa.gov/documentation/approximations/
# Use Higher Order model
md = setflowequation(md, "HO", "all")  # ~10min to run for 500m grid
# Use Full Stokes flow model
md = setflowequation(md, "FS", "all")  # ~35min to run for 500m grid
assert md.flowequation.isFS is True

# Solve
md.toolkits = toolkits()
md.cluster = generic(
    "name", os.uname()[1], "np", 42, "executionpath", os.path.abspath("Models")
)
md = solve(md, "Stressbalance")

# Update model friction fields accordingly
md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient

plotmodel(md, "data", md.friction.coefficient, "data", md.materials.rheology_B)
# Save model
# savevars("Models/pig_control_drag.dat", {"md": md})
# savevars("Models/pig_control_drag_higherorder.dat", {"md": md})
export_netCDF(md=md, filename=f"Models/{md.miscellaneous.name}.nc")

# %%
# Step 5 - Export mesh
import exportVTK as exportVTK  # TODO export to vtk format!

# mdFS = loadmodel("Models/pig_control_drag_fullstokes.dat")
mdFS = loadmodel(path="Models/pig_control_drag_fullstokes.nc")

df = pd.DataFrame(
    data={
        "x": mdFS.mesh.x,
        "y": mdFS.mesh.y,
        "z": mdFS.mesh.z,
        "friction": mdFS.friction.coefficient.squeeze(),
        "velocity": mdFS.results.StressbalanceSolution.Vel.squeeze(),
        "isbasal": mdFS.mesh.vertexonbase,
        "issurface": mdFS.mesh.vertexonsurface,
    },
    columns=["x", "y", "z", "friction", "velocity", "isbasal", "issurface"],
)
df.to_csv("Models/mdFS_xyz_vel_friction.csv", index=False, sep=" ")

# %%
# Step 6 - Plot contour plots
# https://fabrizioguerrieri.com/blog/2017/9/7/surface-graphs-with-irregular-dataset
region = "/".join([str(c) for c in (df.x.min(), df.x.max(), df.y.min(), df.y.max())])
# #!gmt contour Models/mdFS_xyz_vel_friction.csv -Wthin -Cjet -I -Jx1:15000 -R{region} -pdf mdHOvel

df = pd.read_csv("Models/mdFS_xyz_vel_friction.csv", sep=" ")

# Load data
triang = matplotlib.tri.Triangulation(x=df.x, y=df.y)
z = df.velocity

# Plot mesh nodes
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.triplot(
    triang,
    c="#D3D3D3",
    marker=".",
    markerfacecolor="#DC143C",
    markeredgecolor="black",
    markersize=10,
)
ax.set_xlabel("X")
ax.set_ylabel("Y")
plt.show()

# Flot 3D tri surface
fig = plt.figure(figsize=(10, 12))
ax = fig.add_subplot(1, 1, 1, projection="3d")
ax.plot_trisurf(triang, z, cmap="jet")
# ax.scatter(x,y,z, marker='.', s=10, c="black", alpha=0.5)
ax.view_init(elev=60, azim=270)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()


# %% [raw]
# # Inverting for ice rheology (parameter B)
# if step == 8:
#     # Load Model
#     md = loadmodel(path="Pig/Models/PIG_Parameterization_py.dat")
#
#     # Control general
#     md.inversion.iscontrol = 1
#     md.inversion.maxsteps = 40
#     md.inversion.maxiter = 40
#     md.inversion.dxmin = 0.1
#     md.inversion.gttol = 1.0e-6
#     md.verbose = verbose("control", True)
#
#     # Cost functions
#     md.inversion.cost_functions = [101, 103, 502]
#     md.inversion.cost_functions_coefficients = np.ones(
#         shape=(md.mesh.numberofvertices, 3)
#     )
#     md.inversion.cost_functions_coefficients[:1] = 1000
#     md.inversion.cost_functions_coefficients[:2] = 1
#     md.inversion.cost_functions_coefficients[:3] = 1.0e-16
#
#     # Controls
#     md.inversion.control_parameters = ["MaterialsRheologyBbar"]
#     # md.inversion.min_parameters = np.expand_dims(a=md.materials.rheology_B, axis=-1)
#     # md.inversion.max_parameters = np.expand_dims(a=md.materials.rheology_B, axis=-1)
#     # pos = np.argwhere(md.mask.groundedice_levelset < 0)
#     # md.inversion.min_parameters[pos] = cuffey(273)
#     # md.inversion.max_parameters[pos] = cuffey(200)
#     md.inversion.min_parameters = cuffey(273) * np.ones(
#         shape=(md.mesh.numberofvertices, 1)
#     )
#     md.inversion.max_parameters = cuffey(200) * np.ones(
#         shape=(md.mesh.numberofvertices, 1)
#     )
#
#     # Additional parameters
#     md.stressbalance.restol = 0.01
#     md.stressbalance.reltol = 0.1
#     md.stressbalance.abstol = np.NaN
#
#     # Solve
#     md.cluster = generic(
#         "name", os.uname()[1], "np", 2, "executionpath", os.path.abspath("Pig/Models")
#     )
#     # mds = md.extract(md.mask.groundedice_levelset < 0)
#     md = solve(md, "Stressbalance")
#
#     # Update model rheology_B accordingly
#     md.materials.rheology_B = md.results.StressbalanceSolution.MaterialsRheologyBbar
#     plotmodel(md, "data", md.materials.rheology_B)
#
#     # Save model
#     savevars("Pig/Models/PIG_Control_B_py.dat", {"md": md})
#
#     import pandas as pd
#
#     df = pd.DataFrame(
#         data={
#             "x": md.mesh.x,
#             "y": md.mesh.y,
#             "mdSSArheologyB": np.squeeze(md.materials.rheology_B),
#         },
#         columns=["x", "y", "mdSSArheologyB"],
#     )
#     df.to_csv("Pig/Results/mdSSA_xy_rheologyB.csv", index=False, sep=" ")


# %% [markdown]
# ## Plot figures

# %%
# Full Stokes model inverted velocity/friction
dfFS = pd.read_csv("Models/mdFS_xyz_vel_friction.csv", sep=" ")
region = "/".join(
    [str(c) for c in (dfFS.x.min(), dfFS.x.max(), dfFS.y.min(), dfFS.y.max())]
)
# z_attribute = "friction"
z_attribute = "velocity"

if z_attribute == "friction":
    df = dfFS.query(expr="isbasal == True")
if z_attribute == "velocity":
    df = dfFS.query(expr="issurface == True")

x = df.x
y = df.y
z = df[z_attribute]

fig = pygmt.Figure()
pygmt.makecpt(
    cmap="hawaii", series=[z.min(), z.max(), (z.max() - z.min()) / 10], reverse=True,
)
fig.basemap(region=region, projection="x1:1000000", frame=True)
fig.contour(
    x=x, y=y, z=z, I=True, C=True, label_placement="l-1750000/-350000/-1400000/0"
)
fig.colorbar(position="JRM", S=True)
# fig.plot(data="Models/triangle.ijk", pen="thinner")
# fig.savefig(fname=f"Models/{z_attribute}.png")
fig.show()

# %%
# Shallow Shelf Approximation model inverted ice rheology parameter B
dfSSA = pd.read_csv("Pig/Results/mdSSA_xy_rheologyB.csv", sep=" ")
region = "/".join(
    [str(c) for c in (dfSSA.x.min(), dfSSA.x.max(), dfSSA.y.min(), dfSSA.y.max())]
)
z_attribute = "mdSSArheologyB"

x = dfSSA.x
y = dfSSA.y
z = dfSSA[z_attribute]

fig = pygmt.Figure()
pygmt.makecpt(
    cmap="batlow",
    series=[z.min(), z.max(), (z.max() - z.min()) / 100],
    # reverse=True,
)
fig.basemap(region=region, projection="x1:1000000", frame=True)
fig.contour(x=x, y=y, z=z, I=True, C=True)
fig.colorbar(position="JRM", S=True)
fig.savefig(fname=f"Models/{z_attribute}.png")
fig.show()

# %%
# grdtrack on friction/velocity/rheologyB
# !gmt triangulate Pig/Results/_xyz_vel_friction.csv -M -i0,1,4 -h1 -R{region} -I250 -GPig/Results/vel.nc > Pig/Results/triangle.ijk
# !gmt triangulate Pig/Results/mdHO_xyz_vel_friction.csv -M -i0,1,3 -h1 -R{region} -I250 -GPig/Results/friction.nc
# !gmt triangulate Pig/Results/mdSSA_xy_rheologyB.csv -M -i0,1,2 -h1 -R{region} -I250 -GPig/Results/rheologyB.nc

# #!gmt triangulate Models/mdFS_xyz_vel_friction.csv -M -i0,1,4 -h1 -R{region} -I250 -GModels/velocity.nc

z_name = "velocity"  # "rheologyB"
grid = xr.open_dataarray(f"Models/{z_name}.nc")

points = pd.DataFrame(
    data=np.linspace(start=(-1631_500, -259_000), stop=(-1536_500, -95_000), num=50),
    columns=["x", "y"],
)
transect = pygmt.grdtrack(points=points, grid=grid, newcolname=z_name)

transect.plot(x="x", y=z_name)

fig = pygmt.Figure()
pygmt.makecpt(
    cmap="hawaii",
    series=[z.min(), z.max(), (z.max() - z.min()) / 10],
    reverse=True,
    continuous=True,
)
fig.grdimage(grid=grid, cmap=True, frame=True)
fig.plot(x=transect.x, y=transect.y, color=transect[z_name])
fig.colorbar(position="JRM", S=True)
fig.show()
