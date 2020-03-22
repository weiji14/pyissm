import os
import sys

sys.version
sys.path.append("/opt/issm/trunk/bin")
sys.path.append("/opt/issm/trunk/lib")
sys.path.append("/opt/issm/trunk/src/m/dev")

os.chdir("../issm")
os.getcwd()
# sys.path.append("../issm/issm_python/")

assert os.environ["ISSM_DIR"] == "/opt/issm/trunk"
os.environ["PYTHONSTARTUP"] = "/opt/issm/trunk/src/m/dev/devpath.py"

import devpath

# import ISSM

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# %matplotlib inline

from model import model
from bamg import bamg
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

# %%
step = 7

assert os.path.basename(os.getcwd()) == "issm"

# Mesh Generation #1
if step == 1:
    # md = model()
    # Mesh parameters
    domain = "Pig/DomainOutline.exp"
    assert os.path.exists(domain)
    hinit = 10000  # element size for the initial mesh
    hmax = 40000  # maximum element size of the final mesh
    hmin = 5000  # minimum element size of the final mesh
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
    savevars("Pig/Models/PIG_Mesh_generation_py.dat", {"md": md})

# Masks #2
if step == 2:
    md = loadmodel(path="Pig/Models/PIG_Mesh_generation_py.dat")

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
    )  # ice is present when negative

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
    savevars("Pig/Models/PIG_SetMask_py.dat", {"md": md})

# Parameterization #3
if step == 3:
    md = loadmodel(path="Pig/Models/PIG_SetMask_py.dat")

    md = parameterize(md=md, parametername="Pig/Pig_py.par")

    # Use a MacAyeal flow model
    md = setflowequation(md, "SSA", "all")

    # Save model
    savevars("Pig/Models/PIG_Parameterization_py.dat", {"md": md})

# Control Method #4
if step == 4:
    md = loadmodel(path="Pig/Models/PIG_Parameterization_py.dat")

    # Control general
    md.inversion.iscontrol = 1
    md.inversion.maxsteps = 20
    md.inversion.maxiter = 40
    md.inversion.dxmin = 0.1
    md.inversion.gttol = 1.0e-4
    md.verbose = verbose("control", True)

    # Cost functions
    md.inversion.cost_functions = [101, 103, 501]
    md.inversion.cost_functions_coefficients = np.ones(
        shape=(md.mesh.numberofvertices, 3)
    )
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

    # Solve
    md.toolkits = toolkits()
    md.cluster = generic(
        "name", os.uname()[1], "np", 2, "executionpath", os.path.abspath("Pig/Models")
    )
    md = solve(md, "Stressbalance")
    md.inversion.control_parameters

    md.materials.rheology_B
    # Update model friction fields accordingly
    md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient

    plotmodel(md, "data", md.friction.coefficient, "data", md.materials.rheology_B)

    # Save model
    savevars("Pig/Models/PIG_Control_drag_py.dat", {"md": md})

# Plot #5
if step == 5:
    md = loadmodel(path="Pig/Models/PIG_Control_drag_py.dat")

    plotmodel(
        md,
        "data",
        md.initialization.vel,
        "title",
        "Observed velocity",
        "data",
        md.results.StressbalanceSolution.Vel,
        "title",
        "Modeled Velocity",
        "data",
        md.geometry.base,
        "title",
        "Bed elevation",
        "data",
        md.results.StressbalanceSolution.FrictionCoefficient,
        "title",
        "Friction Coefficient",
        "colorbar#all",
        "on",
        "colorbartitle#1-2",
        "(m/yr)",
        "caxis#1-2",
        ([1.5, 4000]),
        "colorbartitle#3",
        "(m)",
        "log#1-2",
        10,
    )

# Higher-Order #6
if step == 6:
    # Load Model
    md = loadmodel(path="Pig/Models/PIG_Control_drag_py.dat")

    # Disable inversion
    md.inversion.iscontrol = 0

    # Extrude Mesh
    md = model.extrude(md, 3, 1)
    assert md.mesh.numberoflayers == 3

    # Set Flowequation
    md = setflowequation(md, "HO", "all")

    # Solve
    md = solve(md, "Stressbalance")

    # Save Model
    savevars("Pig/Models/PIG_ModelHO_py.dat", {"md": md})

# Plot #7
if step == 7:
    mdHO = loadmodel("Pig/Models/PIG_ModelHO_py.dat")
    mdSSA = loadmodel("Pig/Models/PIG_Control_drag_py.dat")

    basal = np.argwhere(mdHO.mesh.vertexonbase)
    surf = np.argwhere(mdHO.mesh.vertexonsurface)

    # procdata = ho_minus_obs_vel
    # datasize = (np.shape(procdata)[0], 1)
    # datasize[1] > 1 and datasize[0] != md.mesh.numberofvertices + 1

    # https://fabrizioguerrieri.com/blog/2017/9/7/surface-graphs-with-irregular-dataset
    import matplotlib

    triang = matplotlib.tri.Triangulation(x=mdHO.mesh.x, y=mdHO.mesh.y)
    z = mdHO.initialization.vel

    np.testing.assert_equal(mdHO.mesh.x, mdSSA.mesh.x)

    if 1 == 1:
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
        # plt.show()

    if 2 == 2:
        fig = plt.figure(figsize=(10, 12))
        ax = fig.add_subplot(1, 1, 1, projection="3d")

        ax.plot_trisurf(triang, z, cmap="jet")
        # ax.scatter(x,y,z, marker='.', s=10, c="black", alpha=0.5)
        ax.view_init(elev=60, azim=270)

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.show()
    mdHO.mesh.vertexonbase
    mdHO.friction
    mdHO.mesh

    import pandas as pd

    df = pd.DataFrame(
        data={
            "x": mdHO.mesh.x,
            "y": mdHO.mesh.y,
            "z": mdHO.mesh.z,
            "mdHOvel": mdHO.initialization.vel,
            "mdHOfriction": mdHO.friction.coefficient.squeeze(),
            "isbasal": mdHO.mesh.vertexonbase,
            "issurface": mdHO.mesh.vertexonsurface,
        },
        columns=["x", "y", "z", "mdHOfriction", "mdHOvel", "isbasal", "issurface"],
    )
    df.to_csv("Pig/Results/mdHO_xyz_vel_friction.csv", index=False, sep=" ")
    region = "/".join(
        [str(c) for c in (df.x.min(), df.x.max(), df.y.min(), df.y.max())]
    )

    #!gmt pscontour Pig/Results/xy_mdHOvel.csv -Wthin -Cjet -I -Jx1:15000 -R{region} > Pig/Results/mdHOvel.ps

    plotmodel(
        mdHO,
        "nrows",
        3,
        "ncols",
        2,
        "axis#all",
        "equal",
        "data",
        mdHO.initialization.vel,
        "title",
        "Observed velocity",
        # "data",
        # (
        #     np.squeeze(mdHO.results.StressbalanceSolution.Vel[surf])
        #     - np.squeeze(mdHO.initialization.vel[surf])
        # ),
        # "title",
        # "(HO-observed) velocities",
        "data",
        mdSSA.results.StressbalanceSolution.Vel,
        "title",
        "Modeled SSA Velocity",
        # "data",
        # (
        #     mdHO.results.StressbalanceSolution.Vel[surf][:, :, 0]
        #     - mdSSA.results.StressbalanceSolution.Vel
        # ),
        # "title",
        # "(HO-SSA) velocities",
        "data",
        mdHO.results.StressbalanceSolution.Vel,
        "title",
        "Modeled HO surface Velocities",
        # "data",
        # (
        #     mdHO.results.StressbalanceSolution.Vel[surf][:, :, 0]
        #     - mdHO.results.StressbalanceSolution.Vel[basal]
        # ),
        # "title",
        # "(HOsurf-HO base) velocities",
        "caxis#1",
        ([1.5, 4000]),
        "caxis#3",
        ([1.5, 4000]),
        "caxis#5",
        ([1.5, 4000]),
        "colorbar#all",
        "on",
        "view#all",
        2,
        "colorbartitle#all",
        "(m/yr)",
        "layer#5",
        1,
        "log#1",
        10,
        "log#3",
        10,
        "log#5",
        10,
    )
