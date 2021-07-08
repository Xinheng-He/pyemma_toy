# pyemma_toy
A toy model for pyemma 2.5.7


Use Python 3.7.3 with the following packages

mdtraj 1.9.3
pyemma 2.5.7
numpy 1.17.3
matplotlib 3.1.2
pandas 0.25.3 

It has been tested on Python 3.7.3, any capable platform supports this demo.

The anticipated outputs are the following：
the position of microstates on the free energy landscape (cluster.png)
implied timescale test (timescale_txt.png),
c-k test (ck-test-bys-3.png), 
proportion of each state (3states-macrostate_proportion.txt),
the position of each macrostate (3-multiple.png),
the transition time between states (3-data.csv),
the representative structures (3-macro-1-0.1.pdb and so on),
and the representative trajectories (3-macro-1-0.1.nc and so on).

To run it, just type python pyemma_hxh.py in your cmd. It needs around 10 minutes to run on a normal computer.

In my simulation data, I read all trajectories to extract structure and do tICA analysis.
