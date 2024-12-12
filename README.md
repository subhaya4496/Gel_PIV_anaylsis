After the data is obtained from the PIV, we process the data to obtain strain heatmaps and plots for radial profiles.
1. strain_calc_trials.m - This is the main file to obtain radial and azimuthal component of the displacement and the strain components from the PIV file with corresponding mask files.
2. strain_eigen.m - The eigenvalues and eigenvectors of the strain matrix is calculated here, from the data that is obtained from the 'strain_calc_trials.m'.
3. strain_u_r.m - The radial and azimuthal displacements and strains are obtained with respect to the center of contraction from the data that is obtained from the 'strain_calc_trials.m'.
4. Calc_urvr.m - The radial profile of the radial displacements is obtained for each frame and plotted here.
5. Calc_strain.m - The radial profiles of the strain is obtained for each frame and plotted here.
6. strain_plots.m - The heatmaps of different components of the strain can be plotted here.
