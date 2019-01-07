### In working directory:
    scons -c && scons --config=force

Stay in working directory. Open another terminal from the terminal in use, execute in the first terminal:

    ./runSolid

And in the second:

    ./FluidSolver/build/NumSim -geom FluidSolver/geom/default.geom -param FluidSolver/param/default.param
