### After git cloning execute in working directory:

(First, on the IPVS cluster execute

    source /usr/local.nfs/sgs/software/NumSim/bashrc.sh

This is not necessary if your preCICE adapter is already executable globally.)

Then, execute

    scons -c && scons --config=force

Stay in working directory.

#### a) Forced convection over a heated plate:

Open another terminal (now there are two terminals open in parallel in working directory). Then execute in terminal 1

    ./FluidSolver/build/NumSim -geom FluidSolver/geom/plate.geom -param FluidSolver/param/plate.param

and in terminal 2

    laplacianFoam -case Solid_plate/

#### b) Natural convection in cavity with heat-conducting walls:

Open another terminal (now there are two terminals open in parallel in working directory). Then execute in terminal 1

    ./FluidSolver/build/NumSim -geom FluidSolver/geom/caviConv.geom -param FluidSolver/param/caviConv.param

and in terminal 2

    laplacianFoam -case Solid_convection/
