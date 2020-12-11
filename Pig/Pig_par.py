import numpy as np
import xarray as xr

from InterpFromGridToMesh import InterpFromGridToMesh
from paterson import paterson
from m1qn3inversion import m1qn3inversion
from SetMarineIceSheetBC import SetMarineIceSheetBC


# Parameters to change/Try
friction_coefficient = 100  # default [10]
Temp_change = 0  # default [0 K]

# Name and Coordinate system
md.miscellaneous.name = "PIG"
md.mesh.epsg = 3031

# NetCdf Loading
print("   Loading DeepBedMap, BedMachine and ALBMAP data from NetCDF")

# with xr.open_dataset("Data/deepbedmap3_big_int16.nc") as ncdata0:
#    x0 = ncdata0.x.data
#    y0 = ncdata0.y.data
#    bed = ncdata0.z.data.astype(np.float64)

with xr.open_dataset("Data/BedMachineAntarctica_2020-07-15_v02.nc") as ncdata1:
    x1 = ncdata1.x.data.astype(np.int64)
    y1 = np.flipud(ncdata1.y.data).astype(np.int64)
    usrf = np.flipud(ncdata1.surface.data)
    # topg = np.flipud(ncdata1.bed.data)

with xr.open_dataset("Data/ALBMAPv1.nc") as ncdata2:
    x2 = ncdata2.x1.data.astype(np.float64)
    y2 = ncdata2.y1.data.astype(np.float64)
    temp = ncdata2.temp.data.astype(np.float64)
    smb = ncdata2.acca.data.astype(np.float64)
    gflux = ncdata2.ghffm.data.astype(np.float64)

# Geometry
print("   Interpolating surface and ice base")
# md.geometry.base = InterpFromGridToMesh(
#     x1, y1, topg, md.mesh.x, md.mesh.y, 0
# )[0]  # BedMachine
# md.geometry.base = InterpFromGridToMesh(x0, y0, bed, md.mesh.x, md.mesh.y, 0)[
#   0
# ]  # DeepBedMap
md.geometry.surface = InterpFromGridToMesh(x1, y1, usrf, md.mesh.x, md.mesh.y, 0)[0]

print("   Constructing thickness")
md.geometry.thickness = md.geometry.surface - md.geometry.base

# ensure hydrostatic equilibrium on ice shelf:
di = md.materials.rho_ice / md.materials.rho_water

# Get the node numbers of floating nodes
pos = np.argwhere(md.mask.ocean_levelset < 0)

# apply a flotation criterion on the precedingly defined nodes and
# redefine base and thickness accordingly
md.geometry.thickness[pos] = (1 / (1 - di)) * md.geometry.surface[pos]
md.geometry.base[pos] = md.geometry.surface[pos] - md.geometry.thickness[pos]
md.geometry.hydrostatic_ratio = np.ones(shape=md.mesh.numberofvertices)  # For Dakota

# Set min thickness to 1 meter
pos0 = np.argwhere(md.geometry.thickness <= 1)
md.geometry.thickness[pos0] = 1
md.geometry.surface = md.geometry.thickness + md.geometry.base
md.geometry.bed = md.geometry.base.copy()
md.geometry.bed[pos] = md.geometry.base[pos] - 1000

np.argwhere(md.geometry.thickness != md.geometry.surface - md.geometry.base)

# Initialization parameters
print("   Interpolating temperatures")
md.initialization.temperature = (
    InterpFromGridToMesh(x2, y2, temp, md.mesh.x, md.mesh.y, 0)[0]
    + 273.15
    + Temp_change
)

print("   Loading velocities data from NetCDF")

with xr.open_dataset("Data/antarctic_ice_vel_phase_map_v01.nc") as nsidc_vel:
    x3 = nsidc_vel.x.data
    y3 = np.flipud(m=nsidc_vel.y.data)
    velx = np.flipud(m=nsidc_vel.VX.data.astype(np.float64))
    vely = np.flipud(m=nsidc_vel.VY.data.astype(np.float64))

print("   Set observed velocities")
md.initialization.vx = InterpFromGridToMesh(x3, y3, velx, md.mesh.x, md.mesh.y, 0)[0]
md.initialization.vy = InterpFromGridToMesh(x3, y3, vely, md.mesh.x, md.mesh.y, 0)[0]
md.initialization.vz = np.zeros(shape=md.mesh.numberofvertices)
md.initialization.vel = np.sqrt(md.initialization.vx ** 2 + md.initialization.vy ** 2)

print("   Set Pressure")
md.initialization.pressure = (
    md.materials.rho_ice * md.constants.g * md.geometry.thickness
)

print("   Construct ice rheological properties")
md.materials.rheology_n = 3 * np.ones(shape=md.mesh.numberofelements)
md.materials.rheology_B = paterson(md.initialization.temperature)

# Forcings
print("   Interpolating surface mass balance")
mass_balance = InterpFromGridToMesh(x2, y2, smb, md.mesh.x, md.mesh.y, 0)[0]
md.smb.mass_balance = mass_balance * md.materials.rho_water / md.materials.rho_ice

print("   Set geothermal heat flux")
md.basalforcings.geothermalflux = InterpFromGridToMesh(
    x2, y2, gflux, md.mesh.x, md.mesh.y, 0
)[0]

# Friction and inversion set up
print("   Construct basal friction parameters")
if hasattr(md.friction, "coefficient"):  # Budd type sliding law parameters
    md.friction.coefficient = friction_coefficient * np.ones(
        shape=md.mesh.numberofvertices
    )
    md.friction.p = 5 * np.ones(shape=md.mesh.numberofelements)
    md.friction.q = 5 * np.ones(shape=md.mesh.numberofelements)
elif hasattr(md.friction, "Cmax"):  # Schoof type sliding law parameters
    md.friction.C = friction_coefficient * np.ones(shape=md.mesh.numberofvertices)
    # Cmax is from Brondex 2019 https://doi.org/10.5194/tc-13-177-2019
    md.friction.Cmax = 0.4 * np.ones(shape=md.mesh.numberofvertices)
    md.friction.m = (1 / 3) * np.ones(shape=(md.mesh.numberofelements, 1))


# no friction applied on floating ice
pos = np.argwhere(md.mask.ocean_levelset < 0)
if hasattr(md.friction, "coefficient"):  # Budd/Weertman friction
    md.friction.coefficient[pos] = 0
elif hasattr(md.friction, "Cmax"):  # Schoof friction
    md.friction.C[pos] = 0
md.groundingline.migration = "SubelementMigration"

md.inversion = m1qn3inversion()
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy
md.inversion.vel_obs = md.initialization.vel

print("   Set boundary conditions")
md = SetMarineIceSheetBC(md)
md.basalforcings.floatingice_melting_rate = np.zeros(shape=md.mesh.numberofvertices)
md.basalforcings.groundedice_melting_rate = np.zeros(shape=md.mesh.numberofvertices)
md.thermal.spctemperature = md.initialization.temperature
md.masstransport.spcthickness = np.NaN * np.ones(shape=md.mesh.numberofvertices)
