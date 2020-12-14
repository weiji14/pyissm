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
