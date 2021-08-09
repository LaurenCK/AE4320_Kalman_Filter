Lauren Kaffa
9-8-2021
Kalman Filter Assignment AE4320
----------------------------------------------------------------------

Code consists of 
1) Data generation: 						data_generation.m (saves in Data and Data2 --> data with higher airdata noise level)
2) Kalman Filter:						Kalman.m (saves in Kalman_data and Kalman_Data2 --> data with higher airdata noise level)
3) Kalman plots:						Kalman_plotting_script.m
4) KF sensitivity analysis:					sensitivity_analysis.m
5) Calculation aerodynamic model coefficients			measurement_vector_generation.m (saves in Measurement_vector_data)
6) Parameter estimation:					parameter_estimation_maneuvers.m (saves in Training_data)
7) Parameter plots:						parameter_estimation_maneuvers_plots.m 
8) Parameter analysis:						parameter_analysis.m
9) Alternative model structure					model_structure_analysis.m
10) Runge-kutta scheme						rk4.m (reference: AE4320 course)
11) calc_dF, calc_dH, calc_G, calc_h, def_f and def_h:		Used in Kalman.m to calculate the system as defined in the assignment. def_f and def_h are used to calculate symbolically the Jacobians df/dx and dh/dx 
12) kf_calc_ObsRank.m						Used in Kalman.m to check whether the system is observable (reference: AE4320 course)							


Explanation:

1) Generation of data to be used in the state estimation (Kalman Filter). Bias and noise terms added. Also plots the results for the da3211 maneuver.
2) Executes the IEKF and takes input the STD of the VTAS, AoA and sideslip angle (airdata measurements) and sensitiviy bool.
   When sensitivity == 1 data is saved in Kalman_Data2.
3) Plots all figures in the KF chapters in report.
4) Calculate and plot the change in KF results due to 10x higher variance V, alpha, beta.
5) Generation of the dimensionless aerodynamic force- and moment components.
6) Seperate training and validation data and train model parameters on the training data using OLS regression.
7) Plot all results from the parameter estimation trained on the training data, including results on validation data and residual analysis plots.
8) Calculate parameter covariance matrices and plot them.
9) Construct new aerodynamic model structures and train on training data. Plots results and residuals on validation data.