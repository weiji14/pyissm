import os

import numpy as np
import pandas as pd
import pygmt
import xarray as xr

os.chdir("../issm")

# %% Higher order model inverted velocity/friction
dfHO = pd.read_csv("Pig/Results/mdHO_xyz_vel_friction.csv", sep=" ")
region = "/".join(
    [str(c) for c in (dfHO.x.min(), dfHO.x.max(), dfHO.y.min(), dfHO.y.max())]
)
z_attribute = "mdHOfriction"
# z_attribute = "mdHOvel"

if z_attribute == "mdHOfriction":
    df = dfHO.query(expr="isbasal == True")
if z_attribute == "mdHOvel":
    df = dfHO.query(expr="issurface == True")

x = df.x
y = df.y
z = df[z_attribute]

fig = pygmt.Figure()
pygmt.makecpt(
    cmap="hawaii",
    series=[str(z.min()), str(z.max()), str((z.max() - z.min()) / 100)],
    reverse=True,
)
fig.basemap(region=region, projection="x1:1000000", frame=True)
fig.contour(
    x=x, y=y, z=z, I=True, C=True, label_placement="l-1750000/-350000/-1400000/0"
)
fig.colorbar(position="JRM", S=True)
# fig.plot(data="Pig/Results/triangle.ijk", pen="thinner")
# fig.savefig(fname=f"Pig/Results/{z_attribute}.png")
fig.show()


# %% Shallow Shelf Approximation model inverted ice rheology parameter B
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
    series=[str(z.min()), str(z.max()), str((z.max() - z.min()) / 100)],
    # reverse=True,
)
fig.basemap(region=region, projection="x1:1000000", frame=True)
fig.contour(x=x, y=y, z=z, I=True, C=True)
fig.colorbar(position="JRM", S=True)
fig.savefig(fname=f"Pig/Results/{z_attribute}.png")
fig.show()


# %% grdtrack on friction/velocity/rheologyB
#!gmt triangulate Pig/Results/mdHO_xyz_vel_friction.csv -M -i0,1,4 -h1 -R{region} -I250 -GPig/Results/vel.nc > Pig/Results/triangle.ijk
#!gmt triangulate Pig/Results/mdHO_xyz_vel_friction.csv -M -i0,1,3 -h1 -R{region} -I250 -GPig/Results/friction.nc
#!gmt triangulate Pig/Results/mdSSA_xy_rheologyB.csv -M -i0,1,2 -h1 -R{region} -I250 -GPig/Results/rheologyB.nc
z_name = "rheologyB"
grid = xr.open_dataarray(f"Pig/Results/{z_name}.nc")

points = pd.DataFrame(
    data=np.linspace(start=(-1750_000, -350_000), stop=(-1400_000, 0), num=50),
    columns=["x", "y"],
)
transect = pygmt.grdtrack(points=points, grid=grid, newcolname=z_name)

transect.plot(x="x", y=z_name)

fig = pygmt.Figure()
pygmt.makecpt(
    cmap="hawaii",
    series=[str(z.min()), str(z.max()), str((z.max() - z.min()) / 100)],
    reverse=True,
)
fig.grdimage(grid=grid, cmap=True)
fig.plot(x=transect.x, y=transect.y, color=transect[z_name])
fig.colorbar(position="JRM", S=True)
fig.show()
