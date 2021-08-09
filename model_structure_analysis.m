%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct new aerodynamic model structures and train on training data.
% Include validation data on the trained new model structure. 
% Residuals are calculated. 
% Print and plot all results including residuals on validation data.
% Code is seperated into  1) symmetric and 2) asymmetric part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sym = load("Training_data/symmetric.mat");
asym = load("Training_data/asymmetric.mat");
N_start = 200;     % remove first 2 seconds due to bias estimate not being accurate
b = 13.3250;       % [m]
c = 1.9910;        % [m]
format short

%%%%%%%%%%%%%%%%%%%%%%%
%% Symmetric models %%
%%%%%%%%%%%%%%%%%%%%%%%
de3211 = load("Measurement_vector_data/de3211");

t_training_end = 16.5;
N_training_end = t_training_end/0.01;

Cx_training = de3211.file.file.Cx(N_start+1:N_training_end);
Cz_training = de3211.file.file.Cz(N_start+1:N_training_end);
Cm_training = de3211.file.file.Cm(N_start+1:N_training_end);

Cx_val = de3211.file.file.Cx(N_start+1:end);
Cz_val = de3211.file.file.Cz(N_start+1:end);
Cm_val = de3211.file.file.Cm(N_start+1:end);

alpha_training = de3211.file.file.alpha_KF(N_start+1:N_training_end);
q_training = de3211.file.file.qm(N_start+1:N_training_end);
V_training_sym = de3211.file.file.Vm_KF(N_start+1:N_training_end);
elevator_training = de3211.file.file.de(N_start+1:N_training_end);
% take average of the two power settings
thrust_training = (de3211.file.file.Tc1(N_start+1:N_training_end)+de3211.file.file.Tc2(N_start+1:N_training_end))./2; 
u_training = de3211.file.file.u_KF(N_start+1:N_training_end);
u_dot_training = de3211.file.file.u_dot(N_start+1:N_training_end);
w_training = de3211.file.file.w_KF(N_start+1:N_training_end);

alpha_dot_training = de3211.file.file.alpha_dot(N_start+1:N_training_end);
q_dot_training = de3211.file.file.q_dot(N_start+1:N_training_end);

% NEW MODEL STRUCTURE %
% CX %%
A_Cx_training  = ones(length(Cx_training),9);
A_Cx_training (:,2:9) = [ alpha_training,alpha_training.^2, (q_training.*c)./V_training_sym, elevator_training, thrust_training,alpha_dot_training,u_training,w_training]; 

% CZ %%
A_Cz_training  = ones(length(Cz_training),10);
A_Cz_training (:,2:10) = [ alpha_training, (q_training.*c)./V_training_sym, elevator_training, thrust_training,alpha_dot_training,u_training,u_dot_training,w_training,q_dot_training]; 

% Cm %%
A_Cm_training  = ones(length(Cm_training),9);
A_Cm_training (:,2:9) =  [ alpha_training, (q_training.*c)./V_training_sym, elevator_training, thrust_training,alpha_dot_training,u_training,w_training,q_dot_training]; 

% parameter estimation % 
theta_hat_Cx = (A_Cx_training' * A_Cx_training)^-1 *A_Cx_training'*Cx_training; 
theta_hat_Cz = (A_Cz_training' * A_Cz_training)^-1 *A_Cz_training'*Cz_training; 
theta_hat_Cm = (A_Cm_training' * A_Cm_training)^-1 *A_Cm_training'*Cm_training;


% Estimate on validation data %
% Cx %%
A_Cx_val = ones(length(de3211.file.file.Cx)-N_start,9);
A_Cx_val(:,2:9) = [ de3211.file.file.alpha_KF(N_start+1:end-1),de3211.file.file.alpha_KF(N_start+1:end-1).^2, (de3211.file.file.qm(N_start+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+1:end-1), de3211.file.file.de(N_start+1:end-1), (de3211.file.file.Tc1(N_start+1:end-1)+de3211.file.file.Tc2(N_start+1:end-1))./2,de3211.file.file.alpha_dot(N_start+1:end),de3211.file.file.u_KF(N_start+1:end-1),de3211.file.file.w_KF(N_start+1:end-1)]; 
Cx_est_val_new = A_Cx_val*theta_hat_Cx;

% CZ %%
A_Cz_val = ones(length(de3211.file.file.Cz)-N_start,10);
A_Cz_val(:,2:10) = [ de3211.file.file.alpha_KF(N_start+1:end-1), (de3211.file.file.qm(N_start+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+1:end-1), de3211.file.file.de(N_start+1:end-1), (de3211.file.file.Tc1(N_start+1:end-1)+de3211.file.file.Tc2(N_start+1:end-1))./2,de3211.file.file.alpha_dot(N_start+1:end),de3211.file.file.u_KF(N_start+1:end-1),de3211.file.file.u_dot(N_start+1:end),de3211.file.file.w_KF(N_start+1:end-1),de3211.file.file.q_dot(N_start+1:end)];
Cz_est_val_new = A_Cz_val*theta_hat_Cz;

% Cm %%
A_Cm_val = ones(length(de3211.file.file.Cm)-N_start,9);
A_Cm_val(:,2:9) = [ de3211.file.file.alpha_KF(N_start+1:end-1), (de3211.file.file.qm(N_start+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+1:end-1), de3211.file.file.de(N_start+1:end-1), (de3211.file.file.Tc1(N_start+1:end-1)+de3211.file.file.Tc2(N_start+1:end-1))./2,de3211.file.file.alpha_dot(N_start+1:end),de3211.file.file.u_KF(N_start+1:end-1),de3211.file.file.w_KF(N_start+1:end-1),de3211.file.file.q_dot(N_start+1:end)];
Cm_est_val_new = A_Cm_val*theta_hat_Cm;

total_Cx_res = sum(abs(sym.file.Cx_residuals))
total_Cz_res = sum(abs(sym.file.Cz_residuals))
total_Cm_res = sum(abs(sym.file.Cm_residuals))

Cx_res_val_new = Cx_val-Cx_est_val_new;
Cz_res_val_new = Cz_val-Cz_est_val_new;
Cm_res_val_new = Cm_val-Cm_est_val_new;

total_Cx_res_new = sum(abs(Cx_res_val_new))
total_Cz_res_new = sum(abs(Cz_res_val_new))
total_Cm_res_new = sum(abs(Cm_res_val_new))

perc1 = (total_Cx_res_new-total_Cx_res)/total_Cx_res*100
perc2 = (total_Cz_res_new-total_Cz_res)/total_Cz_res*100
perc3 = (total_Cm_res_new-total_Cm_res)/total_Cm_res*100


%%%%%%%%%%%%%%%%%%
%% Plotting %%%%%%
%%%%%%%%%%%%%%%%%%
close("all")
t_sym = de3211.file.file.t(N_start+1:end-1);

figure(1)
sgtitle("Model Estimates and Residuals on Validation Data with New Model Structures");
subplot(2,3,1)
plot(t_sym, Cx_val)
hold on
plot(t_sym, sym.file.Cx_est_val)
hold on
plot(t_sym, Cx_est_val_new)
title("C_{X} estimate")
xlabel("Time [s]")
ylabel("[-]")

subplot(2,3,2)
plot(t_sym, Cz_val)
hold on
plot(t_sym, sym.file.Cz_est_val)
hold on
plot(t_sym, Cz_est_val_new)
hold on
title("C_{Z} estimate")
xlabel("Time [s]")
ylabel("[-]")

subplot(2,3,3)
plot(t_sym, Cm_val)
hold on
plot(t_sym, sym.file.Cm_est_val)
hold on
plot(t_sym, Cm_est_val_new)
legend("Measured","Original structure estimate","New structure estimate")
title("C_{m} estimate")
xlabel("Time [s]")
ylabel("[-]")

subplot(2,3,4)
plot(t_sym, sym.file.Cx_residuals)
hold on
plot(t_sym, Cx_res_val_new)
title("C_{X} residuals")
xlabel("Time [s]")
ylabel("[-]")

subplot(2,3,5)
plot(t_sym, sym.file.Cz_residuals)
hold on
plot(t_sym, Cz_res_val_new)
title("C_{Z} residuals")
xlabel("Time [s]")
ylabel("[-]")

subplot(2,3,6)
plot(t_sym, sym.file.Cm_residuals)
hold on
plot(t_sym, Cm_res_val_new)
title("C_{m} residuals")
xlabel("Time [s]")
ylabel("[-]")
legend("Original structure residuals","New structure residuals")


%%%%%%%%%%%%%%%%%%%%%%%
%% Asymmetric models %%
%%%%%%%%%%%%%%%%%%%%%%%

dr3211 = load("Measurement_vector_data/dr3211");
da3211 = load("Measurement_vector_data/da3211");
drdoublet = load("Measurement_vector_data/drdoublet");
dadoublet = load("Measurement_vector_data/dadoublet");

% Concatenate the data of the maneuvers together
% Note that the first 2 seconds fro both maneuvers are discarded due to initial bad bias estimate 
Cy_training = [dr3211.file.file.Cy(N_start+1:end);da3211.file.file.Cy(N_start+1:end)];
Cl_training = [dr3211.file.file.Cl(N_start+1:end);da3211.file.file.Cl(N_start+1:end)];
Cn_training = [dr3211.file.file.Cn(N_start+1:end);da3211.file.file.Cn(N_start+1:end)];

% Add validation data (with drdoublet and dadoublet together) %
Cy_val = [drdoublet.file.file.Cy(N_start+1:end);dadoublet.file.file.Cy(N_start+1:end)];
Cl_val = [drdoublet.file.file.Cl(N_start+1:end);dadoublet.file.file.Cl(N_start+1:end)];
Cn_val = [drdoublet.file.file.Cn(N_start+1:end);dadoublet.file.file.Cn(N_start+1:end)];

beta_training_as = [dr3211.file.file.beta_KF(N_start+1:end-1);da3211.file.file.beta_KF(N_start+1:end-1)];
beta_dot_training_as = [dr3211.file.file.beta_dot(N_start+1:end);da3211.file.file.beta_dot(N_start+1:end)];
p_training_as = [dr3211.file.file.pm(N_start+1:end-1);da3211.file.file.pm(N_start+1:end-1)];
p_dot_training_as = [dr3211.file.file.p_dot(N_start+1:end);da3211.file.file.p_dot(N_start+1:end)];
V_training_as = [dr3211.file.file.Vm_KF(N_start+1:end-1);da3211.file.file.Vm_KF(N_start+1:end-1)];
r_training_as = [dr3211.file.file.rm(N_start+1:end-1);da3211.file.file.rm(N_start+1:end-1)];
r_dot_training_as = [dr3211.file.file.r_dot(N_start+1:end);da3211.file.file.r_dot(N_start+1:end)];
aileron_training_as = [dr3211.file.file.da(N_start+1:end-1);da3211.file.file.da(N_start+1:end-1)];
rudder_training_as = [dr3211.file.file.dr(N_start+1:end-1);da3211.file.file.dr(N_start+1:end-1)];
v_dot_training_as = [dr3211.file.file.v_dot(N_start+1:end);da3211.file.file.v_dot(N_start+1:end)];

% VALIDATION
beta_val = [drdoublet.file.file.beta_KF(N_start+1:end-1);dadoublet.file.file.beta_KF(N_start+1:end-1)];
beta_dot_val = [drdoublet.file.file.beta_dot(N_start+1:end);dadoublet.file.file.beta_dot(N_start+1:end)];
p_val = [drdoublet.file.file.pm(N_start+1:end-1);dadoublet.file.file.pm(N_start+1:end-1)];
p_dot_val = [drdoublet.file.file.p_dot(N_start+1:end);dadoublet.file.file.p_dot(N_start+1:end)];
V_val = [drdoublet.file.file.Vm_KF(N_start+1:end-1);dadoublet.file.file.Vm_KF(N_start+1:end-1)];
r_val = [drdoublet.file.file.rm(N_start+1:end-1);dadoublet.file.file.rm(N_start+1:end-1)];
r_dot_val = [drdoublet.file.file.r_dot(N_start+1:end);dadoublet.file.file.r_dot(N_start+1:end)];
v_dot_val = [drdoublet.file.file.v_dot(N_start+1:end);dadoublet.file.file.v_dot(N_start+1:end)];
aileron_val = [drdoublet.file.file.da(N_start+1:end-1);dadoublet.file.file.da(N_start+1:end-1)];
rudder_val = [drdoublet.file.file.dr(N_start+1:end-1);dadoublet.file.file.dr(N_start+1:end-1)];

%-------------------------------------------------------------------------------------------------------------------%

% NEW STRUCTURE: C_v and Cn_beta_dot added
A_C_y_training = ones(length(Cy_training),10);
A_C_y_training(:,2:10) = [beta_training_as,((p_training_as).*b)./(2*V_training_as),(r_training_as.*b)./(2*V_training_as), aileron_training_as,rudder_training_as,beta_dot_training_as,V_training_as,p_dot_training_as,r_dot_training_as];
A_C_y_val = ones(length(Cy_val),10);
A_C_y_val(:,2:10) = [beta_val,((p_val).*b)./(2*V_val),(r_val.*b)./(2*V_val), aileron_val,rudder_val,beta_dot_val,V_val,p_dot_val,r_dot_val];

A_C_l_training = ones(length(Cl_training),10);
A_C_l_training(:,2:10) = [beta_training_as,((p_training_as).*b)./(2*V_training_as),(r_training_as.*b)./(2*V_training_as), aileron_training_as,rudder_training_as,beta_dot_training_as,V_training_as,p_dot_training_as,r_dot_training_as];
A_C_l_val = ones(length(Cl_val),10);
A_C_l_val(:,2:10) = [beta_val,((p_val).*b)./(2*V_val),(r_val.*b)./(2*V_val), aileron_val,rudder_val,beta_dot_val,V_val,p_dot_val,r_dot_val];

A_C_n_training = ones(length(Cn_training),10);
A_C_n_training(:,2:10) = [beta_training_as,((p_training_as).*b)./(2*V_training_as),(r_training_as.*b)./(2*V_training_as), aileron_training_as,rudder_training_as,beta_dot_training_as,V_training_as,p_dot_training_as,r_dot_training_as];
A_C_n_val = ones(length(Cn_val),10);
A_C_n_val(:,2:10) = [beta_val,((p_val).*b)./(2*V_val),(r_val.*b)./(2*V_val), aileron_val,rudder_val,beta_dot_val,V_val,p_dot_val,r_dot_val];


theta_hat_Cy = (A_C_y_training' * A_C_y_training)^-1 *A_C_y_training'*Cy_training; 
theta_hat_Cl = (A_C_l_training' * A_C_l_training)^-1 *A_C_l_training'*Cl_training; 
theta_hat_Cn = (A_C_n_training' * A_C_n_training)^-1 *A_C_n_training'*Cn_training;
Cy_est_new = A_C_y_val*theta_hat_Cy;
Cl_est_new = A_C_l_val*theta_hat_Cl;
Cn_est_new = A_C_n_val*theta_hat_Cn;

Cy_res_val_new = Cy_val - (Cy_est_new);
Cl_res_val_new  = Cl_val - (Cl_est_new);
Cn_res_val_new  = Cn_val - (Cn_est_new);


total_Cy_res = sum(abs(asym.file.Cy_res_val))
total_Cl_res = sum(abs(asym.file.Cl_res_val))
total_Cn_res = sum(abs(asym.file.Cn_res_val))

total_Cy_res_new = sum(abs(Cy_res_val_new))
total_Cl_res_new = sum(abs(Cl_res_val_new))
total_Cn_res_new = sum(abs(Cn_res_val_new))

perc4 = (total_Cy_res_new-total_Cy_res)/total_Cy_res*100
perc5 = (total_Cl_res_new-total_Cl_res)/total_Cl_res*100
perc6 = (total_Cn_res_new-total_Cn_res)/total_Cn_res*100


%%%%%%%%%%%%%%%%%%
%% Plotting %%%%%%
%%%%%%%%%%%%%%%%%%
t_asym = (0:0.01:(length(Cy_training)-1)*0.01)';

figure(2)
sgtitle("Model Estimates and Residuals on Validation Data with New Model Structures");
subplot(2,3,1)
plot(t_asym, Cy_val)
hold on
plot(t_asym, asym.file.Cy_est_val)
hold on
plot(t_asym, Cy_est_new)
title("C_{Y} estimate")
xlabel("Time [s]")
ylabel("[-]")


subplot(2,3,2)
plot(t_asym, Cl_val)
hold on
plot(t_asym, asym.file.Cl_est_val)
hold on
plot(t_asym, Cl_est_new)
hold on
title("C_{l} estimate")
xlabel("Time [s]")
ylabel("[-]")

subplot(2,3,3)
plot(t_asym, Cn_val)
hold on
plot(t_asym, asym.file.Cn_est_val)
hold on
plot(t_asym, Cn_est_new)
title("C_{n} estimate")
xlabel("Time [s]")
ylabel("[-]")
legend("Measured","Original structure estimate","New structure estimate")


subplot(2,3,4)
plot(t_asym, asym.file.Cy_res_val)
hold on
plot(t_asym, Cy_res_val_new)
title("C_{Y} residuals")
xlabel("Time [s]")
ylabel("[-]")

subplot(2,3,5)
plot(t_asym, asym.file.Cl_res_val)
hold on
plot(t_asym, Cl_res_val_new)
title("C_{l} residuals")
xlabel("Time [s]")
ylabel("[-]")

subplot(2,3,6)
plot(t_asym, asym.file.Cn_res_val)
hold on
plot(t_asym, Cn_res_val_new)
title("C_{n} residuals")
xlabel("Time [s]")
ylabel("[-]")
legend("Original structure residuals","New structure residuals")


%% Generate Table %%
theta_hat_final = [string(theta_hat_Cx);string(theta_hat_Cz);string(theta_hat_Cm);string(theta_hat_Cy);string(theta_hat_Cl);string(theta_hat_Cn)];
parameter_names = ["Cx0","Cxalpha","Cxalpha^2","Cxq","Cxdelta_e","CxTc","Cxalpha_dot","Cxu","Cxw","Cz0","Czalpha","Czq","Czdelta_e","CzTc","Calpha_dot","Czu","Czu_dot","Czw","Czq_dot","Cm0","Cmalpha","Cmq","Cmdelta_e","CmTc","Cmalpha_dot","Cmu","Cmw","Cmq_dot","Cy0","Cybeta","Cyp","Cyr","Cydelta_a","Cy_delta_r","Cybeta_dot","Cyv","Cyp_dot","Cyr_dot","Cl0","Clbeta","Clp","Clr","Cldelta_a","Cl_delta_r","Clbeta_dot","Clv","Clp_dot","Clr_dot","Cn0","Cnbeta","Cnp","Cnr","Cndelta_a","Cn_delta_r","Cnbeta_dot","Cnv","Cnp_dot","Cnr_dot"];
table_theta_hats = table(parameter_names',theta_hat_final);
display(table_theta_hats)


