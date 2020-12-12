# %%
import collections
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
from frictionschoof import frictionschoof
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
hmax = 20000  # maximum element size of the final mesh (originally 40000)
hmin = 250  # minimum element size of the final mesh (originally 5000) TODO 250
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
    md, "hmax", hmax, "hmin", hmin, "gradation", gradation, "field", vel_obs, "err", err
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
md.mask.ocean_levelset = groundedice  # ice is grounded for mask equal one
md.mask.ice_levelset = -1 * np.ones(
    shape=md.mesh.numberofvertices
)  # ice is present when negative

# plotting
plotmodel(
    md,
    "data",
    md.mask.ocean_levelset,
    "title",
    "grounded/floating",
    "data",
    md.mask.ice_levelset,
    "title",
    "ice/no-ice",
)

# Save model
export_netCDF(md=md, filename="Models/pig_setmask.nc")

# %%
# Step 3 - Parameterization
for bedname in ("DeepBedMap", "BedMachine"):
    md = loadmodel(path="Models/pig_setmask.nc")

    # DeepBedMap v1.1 https://doi.org/10.5281/zenodo.4054246
    if bedname == "DeepBedMap":
        with xr.open_dataset("Data/deepbedmap_dem.nc") as ncdata0:
            x = ncdata0.x.data
            y = ncdata0.y.data
            bed = ncdata0.z.data.astype(np.float64)
    # BedMachine v2 https://doi.org/10.5067/E1QL9HFQ7A8M
    if bedname == "BedMachine":
        with xr.open_dataset("Data/BedMachineAntarctica_2020-07-15_v02.nc") as ncdata1:
            x = ncdata1.x.data.astype(np.int64)
            y = np.flipud(ncdata1.y.data).astype(np.int64)
            # usrf = np.flipud(ncdata1.surface.data)
            bed = np.flipud(ncdata1.bed.data)

    md.geometry.base = InterpFromGridToMesh(x, y, bed, md.mesh.x, md.mesh.y, 0)[0]
    md.friction = frictionschoof()  # Set to use Schoof (2005) type sliding law
    md = parameterize(md=md, parametername="Pig/Pig_par.py")

    # Save model
    filepath = f"Models/pig_parameterization_{bedname.lower()}.nc"
    if os.path.exists(path=filepath):
        os.remove(path=filepath)
    export_netCDF(md=md, filename=filepath)


# %%
md.friction

# %%
# Step 4 - Control Method inverting for basal friction
for bedname in ("DeepBedMap", "BedMachine"):
    md = loadmodel(path=f"Models/pig_parameterization_{bedname.lower()}.nc")
    md.miscellaneous.name = f"pig_{bedname.lower()}_control_drag_fullstokes"

    # Extrude Mesh
    md = model.extrude(md, 10, 3)  # 10 layers, with extrusion exponent of 3
    assert md.mesh.numberoflayers == 10
    md.friction.m = md.friction.m[
        ..., np.newaxis
    ]  # shape like (123, 1) instead of (123,)
    assert md.friction.m.shape[-1] == 1

    ## Control general
    md.inversion.iscontrol = 1
    # M1QN3 optimizer parameters
    md.inversion.maxsteps = 30  # maximum number of iterations (gradient computation)
    md.inversion.maxiter = 40  # maximum number of Function evaluation
    md.inversion.dxmin = 0.1  #  convergence criterion: two points less than dxmin from each other (sup-norm) are considered identical
    md.inversion.gttol = 1.0e-4  # gradient relative convergence criterion 2
    md.verbose = verbose("control", True)
    # Toolkit for Advanced Optimization (TAO) parameters TODO?

    # Cost functions
    # Following Morlighem et al., 2013 https://doi.org/10.1002/jgrf.20125
    # Alternatively, follow Kyrke-Smith 2018 at https://doi.org/10.3389/feart.2018.00033
    md.inversion.cost_functions = [
        101,  # Surface Absolute Velocity Misfit
        103,  # Surface Log Velocity Misfit
        501,  # Drag Coefficient Absolute Gradient
    ]
    md.inversion.cost_functions_coefficients = np.ones(
        shape=(md.mesh.numberofvertices, 3)
    )
    md.inversion.cost_functions_coefficients[:, 0] = 1
    md.inversion.cost_functions_coefficients[:, 1] = 100
    md.inversion.cost_functions_coefficients[:, 2] = 1e-7

    # Controls
    md.inversion.control_parameters = ["FrictionC"]  # Friction Coefficient
    md.inversion.min_parameters = 1 * np.ones(shape=(md.mesh.numberofvertices, 1))
    md.inversion.max_parameters = 200 ** 2 * np.ones(
        shape=(md.mesh.numberofvertices, 1)
    )  # max friction coefficient C^2

    # Additional parameters
    md.stressbalance.restol = 0.01
    md.stressbalance.reltol = 0.1
    md.stressbalance.abstol = np.NaN

    # Set Flow Equation
    # https://issm.jpl.nasa.gov/documentation/approximations/
    # TODO Use Full Stokes flow model
    md = setflowequation(md, "L1L2", "all")
    # assert md.flowequation.isFS is True

    # Solve
    md.toolkits = toolkits()
    md.cluster = generic(
        "name", os.uname()[1], "np", 4, "executionpath", os.path.abspath("Models")
    )
    md = solve(md, "Stressbalance")

    # Update model friction fields accordingly
    md.friction.coefficient = md.results.StressbalanceSolution.FrictionC
    md.friction.effective_pressure = md.results.StressbalanceSolution.Pressure
    J_mis, J_reg1, J_reg2, J = md.results.StressbalanceSolution.J[-1]

    plotmodel(
        md,
        "data",
        md.friction.coefficient,
        "data",
        md.materials.rheology_B,
        "data",
        md.friction.effective_pressure,
    )
    # Save model
    filepath = f"Models/{md.miscellaneous.name}.nc"
    if os.path.exists(path=filepath):
        os.remove(path=filepath)
    export_netCDF(md=md, filename=filepath)

# %%
# Export pressure/velocity/friction to CSV table
for bedname in ("DeepBedMap", "BedMachine"):
    mdFS = loadmodel(path=f"Models/pig_{bedname.lower()}_control_drag_fullstokes.nc")
    mdFS.results.StressbalanceSolution

    # ISSM regularized Coulomb/Schoof-type sliding law
    # Note that ISSM's Schoof friction coefficient C is actually C ** 2
    C = np.sqrt(mdFS.friction.coefficient.squeeze())  # friction coefficient
    Cmax = mdFS.friction.Cmax.squeeze()  # Iken's bound
    N = mdFS.friction.effective_pressure.squeeze()  # effective pressure
    u_b = mdFS.results.StressbalanceSolution.Vel.squeeze()  # basal velocity
    m = 1 / 3  # power law exponent
    np.testing.assert_equal(actual=mdFS.friction.m, desired=1 / 3)

    # Calculate basal drag \tau_b
    # Get formula from e.g. Brondex et al. 2019 eq.9 https://doi.org/10.5194/tc-13-177-2019
    tau_b = -(
        (C * (abs(u_b) ** (m - 1)) * u_b)
        / ((1 + ((C / (Cmax * N)) ** (1 / m)) * abs(u_b)) ** m)
    )

    df = pd.DataFrame(
        data={
            "x": mdFS.mesh.x,
            "y": mdFS.mesh.y,
            "z": mdFS.mesh.z,
            "pressure": N,
            "friction": C,
            "basal_drag": tau_b.squeeze(),
            "velocity": mdFS.results.StressbalanceSolution.Vel.squeeze(),
            "isbasal": mdFS.mesh.vertexonbase,
            "issurface": mdFS.mesh.vertexonsurface,
        }
    )
    df.to_csv(
        f"Models/mdFS_{bedname.lower()}_xyz_pressure_vel_friction.csv",
        index=False,
        sep=" ",
    )

# %%
# Step 5 - Inverting for ice rheology (parameter B)
# Load Model
bedname = "DeepBedMap"  # "BedMachine"  #
md = loadmodel(path=f"Models/pig_parameterization_{bedname.lower()}.nc")
md.miscellaneous.name = f"pig_{bedname.lower()}_control_rheology_fullstokes"

# Extrude Mesh
md = model.extrude(md, 3, 3)  # TODO 9 layers, with extrusion exponent of 3
assert md.mesh.numberoflayers == 3

# Control general
md.inversion.iscontrol = 1
md.inversion.maxsteps = 40
md.inversion.maxiter = 40
md.inversion.dxmin = 0.1
md.inversion.gttol = 1.0e-6
md.verbose = verbose("control", True)

# Cost functions
md.inversion.cost_functions = [101, 103, 502]
md.inversion.cost_functions_coefficients = np.ones(shape=(md.mesh.numberofvertices, 3))
md.inversion.cost_functions_coefficients[:1] = 1000
md.inversion.cost_functions_coefficients[:2] = 1
md.inversion.cost_functions_coefficients[:3] = 1.0e-16

# Controls
md.inversion.control_parameters = ["MaterialsRheologyBbar"]
# md.inversion.min_parameters = np.expand_dims(a=md.materials.rheology_B, axis=-1)
# md.inversion.max_parameters = np.expand_dims(a=md.materials.rheology_B, axis=-1)
# pos = np.argwhere(md.mask.groundedice_levelset < 0)
# md.inversion.min_parameters[pos] = cuffey(273)
# md.inversion.max_parameters[pos] = cuffey(200)
md.inversion.min_parameters = cuffey(273) * np.ones(shape=(md.mesh.numberofvertices, 1))
md.inversion.max_parameters = cuffey(200) * np.ones(shape=(md.mesh.numberofvertices, 1))

# Additional parameters
md.stressbalance.restol = 0.01
md.stressbalance.reltol = 0.1
md.stressbalance.abstol = np.NaN

# Set Flow Equation
# https://issm.jpl.nasa.gov/documentation/approximations/
# Use Full Stokes flow model
md = setflowequation(md, "FS", "all")  # ~35min to run for 500m grid
assert md.flowequation.isFS is True

# Solve
md.toolkits = toolkits()
md.cluster = generic(
    "name", os.uname()[1], "np", 42, "executionpath", os.path.abspath("Models")
)
# mds = md.extract(md.mask.groundedice_levelset < 0)
md = solve(md, "Stressbalance")

# Update model rheology_B accordingly
md.materials.rheology_B = md.results.StressbalanceSolution.MaterialsRheologyBbar
plotmodel(md, "data", md.materials.rheology_B)

# Save model
export_netCDF(md=md, filename=f"Models/{md.miscellaneous.name}.nc")

# Export rheology to CSV table
mdFS = loadmodel(path=f"Models/pig_{bedname.lower()}_control_rheology_fullstokes.nc")
df = pd.DataFrame(
    data={
        "x": mdFS.mesh.x,
        "y": mdFS.mesh.y,
        "z": mdFS.mesh.z,
        "rheology": np.squeeze(mdFS.materials.rheology_B),
        "isbasal": mdFS.mesh.vertexonbase,
        "issurface": mdFS.mesh.vertexonsurface,
    },
    columns=["x", "y", "z", "rheology", "isbasal", "issurface"],
)
df.to_csv(f"Models/mdFS_{bedname.lower()}_xyz_rheology.csv", index=False, sep=" ")


# %%
# Step 6 - Export mesh
from exportVTK import exportVTK  # TODO export to vtk format!``

# exportVTK(filename="vtkfiles", md=md)

# %% [markdown]
# ## Plot figures

# %%
# Plotting Full Stokes model inverted velocity/friction/rheology
import collections
import dataclasses
import numpy as np
import pandas as pd
import pygmt
import xarray as xr


@dataclasses.dataclass(frozen=True)
class Z_attr:
    varname: str  # z_variable name
    symbol: str  # mathematical symbol
    unit: str  # SI unit of variable
    isbasal: bool  # True if this is a basal variable, False if it's a surface variable


bedname = "DeepBedMap"  # "BedMachine"  #
# z_attr = Z_attr(varname="velocity", symbol="u@-b@-", unit="m yr@-1@-", isbasal=False)
z_attr = Z_attr(
    varname="slipperiness", symbol="C", unit="(Pa yr/m)@+1/2@+", isbasal=True
)
# z_attr = Z_attr(varname="rheology", symbol="B", unit="", isbasal=True)
# z_attr = Z_attr(varname="basal_drag", symbol="@~t@~@-b@-", unit="Pa", isbasal=True)
# z_attr = Z_attr(varname="pressure", symbol="N", unit="Pa", isbasal=True)


if z_attr.varname == "rheology":
    dfFS = pd.read_csv(f"Models/mdFS_{bedname.lower()}_xyz_rheology.csv", sep=" ")
else:
    dfFS = pd.read_csv(
        f"Models/mdFS_{bedname.lower()}_xyz_pressure_vel_friction.csv", sep=" "
    )
    dfFS = dfFS.rename(columns=dict(friction="slipperiness"))
expr: str = "isbasal == True" if z_attr.isbasal else "issurface == True"
df: pd.DataFrame = dfFS.query(expr=expr)[["x", "y", z_attr.varname]]

# df.plot(x="slipperiness", y="velocity", kind="scatter", loglog=False)
# df.plot(x="pressure", y="velocity", kind="scatter", loglog=False)

# %%
# Contour plots of velocity/slipperiness/rheology
xmin, xmax, ymin, ymax, zmin, zmax = pygmt.info(table=df, per_column=True)
region = "/".join(str(i) for i in [xmin, xmax, ymin, ymax])

fig = pygmt.Figure()
pygmt.makecpt(cmap="hawaii", series=[zmin, zmax, (zmax - zmin) / 10], reverse=True)
fig.basemap(
    region=region,
    projection="x1:1000000",
    frame=["af", f'WSne+t"{bedname} {z_attr.varname} {z_attr.symbol}"'],
)
fig.contour(
    data=df.to_numpy(),
    I=True,
    levels=True,
    label_placement="l-1750000/-350000/-1400000/0",
)
fig.colorbar(position="JRM", frame=["af", f'y+l"{z_attr.unit}"'], S=True, xshift="0.5c")
# fig.plot(data="Models/triangle.ijk", pen="thinner")
# fig.savefig(fname=f"Models/{bedname.lower()}_{z_attr.varname}.png")
fig.show()

# %%
# !gmt triangulate Pig/Results/mdFS_xyz_vel_friction.csv -M -i0,1,4 -h1 -R{region} -I250 -GPig/Results/vel.nc > Pig/Results/triangle.ijk
def triangulate(data, **kwargs):
    """Thin wrapper around https://docs.generic-mapping-tools.org/latest/triangulate.html"""
    kind = pygmt.helpers.data_kind(data=data)
    with pygmt.clib.Session() as lib:
        if kind == "file":
            file_context = dummy_context(data)
        elif kind == "matrix":
            file_context = lib.virtualfile_from_matrix(data.values)
        elif kind == "vectors":
            file_context = lib.virtualfile_from_vectors(x, y, z)

        with file_context as infile:
            arg_str = " ".join([infile, pygmt.helpers.build_arg_string(kwargs)])
            lib.call_module(module="triangulate", args=arg_str)


triangulate(data=df, R=region, I=250, G=f"Models/{bedname.lower()}_{z_attr.varname}.nc")

# %%
# Gridded plots, and transect plot of slipperiness/velocity/rheologyB
grid = xr.open_dataarray(f"Models/{bedname.lower()}_{z_attr.varname}.nc")

pointXY = collections.namedtuple(typename="pointXY", field_names="x y")
pointA = pointXY(x=-1590_000, y=-99_000)
pointB = pointXY(x=-1580_000, y=-255_000)
points = pd.DataFrame(
    data=np.linspace(start=pointA, stop=pointB, num=50), columns=["x", "y"]
)
transect = pygmt.grdtrack(points=points, grid=grid, newcolname=z_attr.varname)

fig = pygmt.Figure()
# Plot 2D grid and transect line on map
pygmt.makecpt(
    cmap="hawaii",
    series=[zmin, zmax, (zmax - zmin) / 10],
    reverse=True,
    continuous=True,
)
fig.grdimage(grid=grid, cmap=True, frame=["af", "WSne"], projection="X8c/12c")
fig.text(
    position="TR",
    text=rf"@!\257{z_attr.symbol}={grid.mean().item():.4e}",
    justify="TR",
    offset="-0.2c/-0.2c",
)
fig.plot(x=transect.x, y=transect.y, color=transect[z_attr.varname], pen="1p")
fig.text(x=pointA.x, y=pointA.y, text="A", justify="TR")
fig.text(x=pointB.x, y=pointB.y, text="B", justify="BL")
fig.colorbar(
    position="JRM",
    S=True,
    frame=["1000a" if z_attr.varname == "basal_drag" else "", f'y+l"{z_attr.unit}"'],
    xshift="0.5c",
)
# Plot transect line graph
fig.plot(
    x=transect.x.values,
    y=transect[z_attr.varname].values,
    region=[
        transect.x.min(),
        transect.x.max(),
        transect[z_attr.varname].min(),
        transect[z_attr.varname].max(),
    ],
    projection="X8c/4c",
    frame=["af", f'WSne+t"{bedname} {z_attr.varname}"'],
    pen="2p",
    Y="13c",
)
fig.savefig(fname=f"Models/{bedname.lower()}_{z_attr.varname}_transect.png")
fig.show()
