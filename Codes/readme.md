The order to reproduce all the results should be:
1. effort_allcells.cc Computes the fishing effort in all the cells
2. filtereez.cc Selects high seas grid cells
3. matrixtrajectories.cc Connections between locations (high seas grid cells and ports) following vessels trajectories
4. Apply Infomap to the output of 3, to obtain the fishing provinces.
5. effortincom.cc Distributes the effort in the fishing provinces among the ports.
