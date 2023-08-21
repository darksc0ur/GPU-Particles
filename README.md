# GPU-Particles
A few short MATLAB scripts which demonstrate how easy it is to speed up MATLAB code using a GPU.

just_run_parts_red_and_psi.m <br>
This program solves for multiple particles positions in a flow, where the flow is defined by a streamfunction. This also includes random motion for the particles.

JulianDist2.m <br>
This program searches for particles that have interacted. This program does this by considering all particles lower than a physical distance to have interacted. By binning particles then searching for particles appearing in the same bin.

low_mem_dist.m <br>
This program searches for particles that have interacted. This program does this by considering all particles lower than a physical distance to have interacted. The particles are first sorted into bins. Then a bin is selected and the physical distance between all the particles in that bin, and all of the surrounding bins is calculated.

low_mem_test_sort.m <br>
This program searches for particles that have interacted. This program does this by considering all particles lower than a physical distance to have interacted. This is done by using a binary search array to sort the particle into bins, then computing the distance between the particles in each bin.
