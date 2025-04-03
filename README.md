# kpdr_oscillations

These are MATLAB script files providing example code to produce numerical solutions that appear in 'Stabilising millennial oscillations in large-scale ocean circulation with a delayed feedback due to a circumpolar current' by Andrew Keane, Alexandre Pohl, Henk A. Dijkstra and Andy Ridgwell. This paper has been submitted to Physica D.

The code requires installation of the continuation package DDE-Biftool
https://github.com/DDE-BifTool/DDE-Biftool

The file ocean_box_rhs provides the model in the form of a set of delay differential equations. The file ocean_box_analyse provides examples of:
- Calculating steady-state solutions in one-parameter space and their stability (i.e. eigenvalues of the Jacobian matrix at the steady states);
- Identifies Hopf bifurcations and calculates branches of periodic orbits emerging from the Hopf bifurcations in one-parameter space;
- Calculates curves of Hopf bifurcations in two-parameter space;
- Compute first Lyapunov coefficients of the Hopf bifurcations to determine criticality.
