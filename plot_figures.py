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
            file_context = lib.virtualfile_from_matrix(data.to_numpy())
        elif kind == "vectors":
            file_context = lib.virtualfile_from_vectors(x, y, z)

        with file_context as infile:
            arg_str = " ".join([infile, pygmt.helpers.build_arg_string(kwargs)])
            lib.call_module(module="triangulate", args=arg_str)


triangulate(
    data=df.astype(np.float32),
    R=region,
    I=250,
    G=f"Models/{bedname.lower()}_{z_attr.varname}.nc",
)

# %%
# DeepBedMap z_grid SUBtracted by BedMachine z_grid
with pygmt.clib.Session() as session:
    args = (
        f"Models/deepbedmap_{z_attr.varname}.nc "
        f"Models/bedmachine_{z_attr.varname}.nc "
        f"SUB = Models/diff_{z_attr.varname}.nc "
    )
    session.call_module(module="grdmath", args=args)

# %%
# Gridded plots, and transect plot of slipperiness/velocity/rheologyB
pointXY = collections.namedtuple(typename="pointXY", field_names="x y")
pointA = pointXY(x=-1590_000, y=-99_000)
pointB = pointXY(x=-1580_000, y=-255_000)
points = pd.DataFrame(
    data=np.linspace(start=pointA, stop=pointB, num=250), columns=["x", "y"]
)

# %%
fig = pygmt.Figure()
# Same colormap for both grids
grids = [
    f"Models/{bedname.lower()}_{z_attr.varname}.nc"
    for bedname in ("DeepBedMap", "BedMachine")
]
pygmt.grdinfo(grid=grids[1], T=10)[2:-1]
pygmt.makecpt(
    cmap="hawaii",
    series=pygmt.grdinfo(grid=" ".join(grids), T=10)[2:-1],
    reverse=True,
    continuous=True,
)
with pygmt.clib.Session() as session:

    subplot = lambda args: session.call_module(module="subplot", args=args)
    # Begin subplot
    subplot(args="begin 2x3 -BWSne -Bxaf -Byaf -Fs10c/4c,14c -M0c/0.5c -SRl")

    # Plot transect line graph
    subplot(args=f"set 0,0")
    transectproj = "X32c/4c"
    _xyz = pd.concat(
        pygmt.grdtrack(
            points=points,
            grid=grid,
            newcolname=z_attr.varname,
        )
        for grid in grids
    )
    fig.basemap(
        region=pygmt.info(
            table=_xyz[["x", z_attr.varname]], per_column=True, spacing=10
        ),
        projection=transectproj,
        frame=[
            "af",
            f'WSne+t"{z_attr.varname.title().replace("_", " ")} ({z_attr.symbol}) at Pine Island Glacier"',
        ],
    )
    for bedname in ("DeepBedMap", "BedMachine"):
        grid = f"Models/{bedname.lower()}_{z_attr.varname}.nc"
        transect = pygmt.grdtrack(points=points, grid=grid, newcolname=z_attr.varname)
        fig.plot(
            x=transect.x.values,
            y=transect[z_attr.varname].values,
            projection=transectproj,
            style="c0.1c",
            color="purple" if bedname == "DeepBedMap" else "green",
            label=bedname,
        )
    fig.legend(position="jMR+jMR+o0.2c", S=2, projection=transectproj)
    fig.text(position="TL", text="A", offset="j0.1c", projection=transectproj)
    fig.text(position="TR", text="B", offset="j0.1c", projection=transectproj)

    # Plot 2D grid and transect line on map
    for column, bedname in enumerate(iterable=("DeepBedMap", "BedMachine", "diff")):
        grid = f"Models/{bedname.lower()}_{z_attr.varname}.nc"
        xrgrid = xr.open_dataarray(grid)

        mean_zval = xrgrid.mean().item()
        mean_zval = f"{mean_zval:.4f}" if abs(mean_zval) < 1000 else f"{mean_zval:.4e}"
        title = rf"{bedname}, @!\257{z_attr.symbol}={mean_zval}"
        subplot(args=f"set 1,{column}")  # -A"{title}"

        if bedname == "diff":
            pygmt.makecpt(
                cmap="vik+h0",
                series=pygmt.grdinfo(grid=grid, T="+a0.5+s")[2:-1],
            )
        fig.grdimage(
            grid=grid,
            region=pygmt.grdinfo(grid, I="r")[2:-1],
            cmap=True,
            projection="X10c/14c",
        )
        if bedname == "diff":
            fig.colorbar(position="JMR+o1.5c/0c+e", frame=["af", f'y+l"{z_attr.unit}"'])
        with pygmt.config(FONT_TITLE="14p"):
            fig.plot(
                x=transect.x,
                y=transect.y,
                color=transect[z_attr.varname],
                pen="1p",
                frame=[f'lrbt+t"{title}"'],
            )
        fig.text(x=pointA.x, y=pointA.y, text="A", justify="TR")
        fig.text(x=pointB.x, y=pointB.y, text="B", justify="BL")

    # End subplot
    subplot(args="end")

fig.colorbar(
    position="JBC",
    S=True,
    frame=[
        # "1000af" if z_attr.varname == "basal_drag" else "",
        "xaf",
        f'y+l"{z_attr.unit}"',
    ],
    xshift="2c",
)
fig.savefig(fname=f"Models/inverted_bed_{z_attr.varname}.png")
fig.show()
