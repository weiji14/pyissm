# PyISSM

Experimental code on running the [Ice Sheet System Model (ISSM)](https://issm.jpl.nasa.gov) via Python instead of the default Matlab interface.
Focusing on Pine Island Glacier (based on this [tutorial](https://issm.jpl.nasa.gov/documentation/tutorials/pig/)).

## TODO

- [x] Use newer dataset inputs (e.g. BedMachine, ALBMAPv1, MEaSUREs Phase Map of Antarctic Ice Velocity)
- [ ] Go from Python 2.7 to Python 3
- [ ] Get reproducible [binder](https://mybinder.readthedocs.io) build to work (using Dockerfile)
- [ ] Jupyter notebooks, [PyGMT](https://pygmt.org) plots, and so on!

## Notes

- ISSM currently installed by sysadmin in /opt/issm, activated by running `need issmpy`
- Pine Island Tutorial located in /opt/issm/trunk/examples/Pig, copied to /home/user/pyissm/Pig
- Datasets found under /home/user/pyissm/Data

## References

- Larour, E., Seroussi, H., Morlighem, M., & Rignot, E. (2012). Continental scale, high order, high spatial resolution, ice sheet modeling using the Ice Sheet System Model (ISSM). Journal of Geophysical Research: Earth Surface, 117(F1). https://doi.org/10.1029/2011JF002140
