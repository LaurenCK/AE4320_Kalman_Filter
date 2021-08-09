%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter estimation with data distribution: 1) training data 2) Validation data %%
% Cx, Cz, Cm:
%    % training data: de3211.mat 2-16.5 seconds
%    % validation data: de3211.mat 16.5 - 120 seconds
% Cy Cl, Cn:
%    % training data: dr3211 and da3211
%    % validation data: drdoublet and dadoublet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N_start = 200;     % remove first 2 seconds due to bias estimate not being accurate
b = 13.3250;       % [m]
c = 1.9910;        % [m]

%%%%%%%%%%%%%%%%
%% Cx, Cz, Cm %%
%%%%%%%%%%%%%%%%
de3211 = load("Measurement_vector_data/de3211");

t_training_end = 16.5;
N_training_end = t_training_end/0.01;

file.Cx_training = de3211.file.file.Cx(N_start+1:N_training_end);
file.Cz_training = de3211.file.file.Cz(N_start+1:N_training_end);
file.Cm_training = de3211.file.file.Cm(N_start+1:N_training_end);

file.Cx_val = de3211.file.file.Cx(N_start+1:end);
file.Cz_val = de3211.file.file.Cz(N_start+1:end);
file.Cm_val = de3211.file.file.Cm(N_start+1:end);

alpha_training = de3211.file.file.alpha_KF(N_start+1:N_training_end);
q_training = de3211.file.file.qm(N_start+1:N_training_end);
V_training_sym = de3211.file.file.Vm_KF(N_start+1:N_training_end);
elevator_training = de3211.file.file.de(N_start+1:N_training_end);
% take average of the two power settings
thrust_training = (de3211.file.file.Tc1(N_start+1:N_training_end)+de3211.file.file.Tc2(N_start+1:N_training_end))./2; 

% Cx %%
% theta_Cx = [CX0,CXa,CXa^2,CXq, CXdelta_e,CXT_C]
file.A_Cx_training  = ones(length(file.Cx_training),6);
file.A_Cx_training (:,2:6) = [ alpha_training,alpha_training.^2, (q_training.*c)./V_training_sym, elevator_training, thrust_training]; 

% CZ %%
%theta_Cz = [CZ0,CZa,Czq,CZdelta_e,CZT_C]
file.A_Cz_training  = ones(length(file.Cz_training),5);
file.A_Cz_training (:,2:5) = [ alpha_training, (q_training.*c)./V_training_sym, elevator_training, thrust_training]; 

% Cm %%
%theta_Cm = [Cm0,Cma,Cmq,Cmdelta_e,CmT_c]
file.A_Cm_training  = ones(length(file.Cm_training),5);
file.A_Cm_training (:,2:5) =  [ alpha_training, (q_training.*c)./V_training_sym, elevator_training, thrust_training]; 

% parameter estimation % 
file.theta_hat_Cx = (file.A_Cx_training' * file.A_Cx_training)^-1 *file.A_Cx_training'*file.Cx_training; 
file.theta_hat_Cz = (file.A_Cz_training' * file.A_Cz_training)^-1 *file.A_Cz_training'*file.Cz_training; 
file.theta_hat_Cm = (file.A_Cm_training' * file.A_Cm_training)^-1 *file.A_Cm_training'*file.Cm_training; 

% Estimate on training data (<16) only %
file.Cx_est_train = file.A_Cx_training*file.theta_hat_Cx;
file.Cz_est_train = file.A_Cz_training*file.theta_hat_Cz;
file.Cm_est_train = file.A_Cm_training*file.theta_hat_Cm;

file.Cx_residuals_train = file.Cx_training - file.Cx_est_train;
file.Cz_residuals_train = file.Cz_training - file.Cz_est_train;
file.Cm_residuals_train = file.Cm_training - file.Cm_est_train;

%% Symmetric Training+VALIDATION Data Cx, Cz, Cm: de3211.mat %% 

% Cx %%
% theta_Cx = [CX0,CXa,CXa^2,CXq, CXdelta_e,CXT_C]
file.A_Cx_val = ones(length(de3211.file.file.Cx)-N_start,6);
file.A_Cx_val(:,2:6) = [ de3211.file.file.alpha_KF(N_start+1:end-1),de3211.file.file.alpha_KF(N_start+1:end-1).^2, (de3211.file.file.qm(N_start+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+1:end-1), de3211.file.file.de(N_start+1:end-1), (de3211.file.file.Tc1(N_start+1:end-1)+de3211.file.file.Tc2(N_start+1:end-1))./2]; 
file.Cx_est_val = file.A_Cx_val*file.theta_hat_Cx;

% CZ %%
%theta_Cz = [CZ0,CZa,Czq,CZdelta_e,CZT_C]
file.A_Cz_val = ones(length(de3211.file.file.Cz)-N_start,5);
file.A_Cz_val(:,2:5) = [ de3211.file.file.alpha_KF(N_start+1:end-1), (de3211.file.file.qm(N_start+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+1:end-1), de3211.file.file.de(N_start+1:end-1), (de3211.file.file.Tc1(N_start+1:end-1)+de3211.file.file.Tc2(N_start+1:end-1))./2];
file.Cz_est_val = file.A_Cz_val*file.theta_hat_Cz;

% Cm %%
%theta_Cm = [Cm0,Cma,Cmq,Cmdelta_e,CmT_c]
file.A_Cm_val = ones(length(de3211.file.file.Cm)-N_start,5);
file.A_Cm_val(:,2:5) = [ de3211.file.file.alpha_KF(N_start+1:end-1), (de3211.file.file.qm(N_start+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+1:end-1), de3211.file.file.de(N_start+1:end-1), (de3211.file.file.Tc1(N_start+1:end-1)+de3211.file.file.Tc2(N_start+1:end-1))./2];
file.Cm_est_val = file.A_Cm_val*file.theta_hat_Cm;

% Also store for validation data only after 16 seconds! %
%

% Cx %%
% theta_Cx = [CX0,CXa,CXa^2,CXq, CXdelta_e,CXT_C]
file.A_Cx_val_16 = ones(length(de3211.file.file.Cx)- N_training_end-N_start,6);
file.A_Cx_val_16(:,2:6) = [ de3211.file.file.alpha_KF(N_start+N_training_end+1:end-1),de3211.file.file.alpha_KF(N_start+N_training_end+1:end-1).^2, (de3211.file.file.qm(N_start+N_training_end+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+N_training_end+1:end-1), de3211.file.file.de(N_start+N_training_end+1:end-1), (de3211.file.file.Tc1(N_start+N_training_end+1:end-1)+de3211.file.file.Tc2(N_start+N_training_end+1:end-1))./2]; 
file.Cx_est_val_16 = file.A_Cx_val_16*file.theta_hat_Cx;
file.Cx_val_16 = file.Cx_val(N_training_end+1:end);
file.Cx_residuals_16 = file.Cx_val_16-file.Cx_est_val_16;

% CZ %%
%theta_Cz = [CZ0,CZa,Czq,CZdelta_e,CZT_C]
file.A_Cz_val_16 = ones(length(de3211.file.file.Cz)-N_start-N_training_end,5);
file.A_Cz_val_16(:,2:5) = [ de3211.file.file.alpha_KF(N_start+N_training_end+1:end-1), (de3211.file.file.qm(N_start+N_training_end+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+N_training_end+1:end-1), de3211.file.file.de(N_start+N_training_end+1:end-1), (de3211.file.file.Tc1(N_start+N_training_end+1:end-1)+de3211.file.file.Tc2(N_start+N_training_end+1:end-1))./2];
file.Cz_est_val_16 = file.A_Cz_val_16*file.theta_hat_Cz;
file.Cz_val_16 = file.Cz_val(N_training_end+1:end);
file.Cz_residuals_16 = file.Cz_val_16-file.Cz_est_val_16;

% Cm %%
%theta_Cm = [Cm0,Cma,Cmq,Cmdelta_e,CmT_c]
file.A_Cm_val_16 = ones(length(de3211.file.file.Cm)-N_start-N_training_end,5);
file.A_Cm_val_16(:,2:5) = [ de3211.file.file.alpha_KF(N_start+N_training_end+1:end-1), (de3211.file.file.qm(N_start+N_training_end+1:end-1).*c)./de3211.file.file.Vm_KF(N_start+N_training_end+1:end-1), de3211.file.file.de(N_start+N_training_end+1:end-1), (de3211.file.file.Tc1(N_start+N_training_end+1:end-1)+de3211.file.file.Tc2(N_start+N_training_end+1:end-1))./2];
file.Cm_est_val_16 = file.A_Cm_val_16*file.theta_hat_Cm;
file.Cm_val_16 = file.Cm_val(N_training_end+1:end);
file.Cm_residuals_16 = file.Cm_val_16- file.Cm_est_val_16;


%% Store Symmetric Residuals %%
file.Cx_residuals = de3211.file.file.Cx(N_start+1:end) - file.Cx_est_val;
file.Cz_residuals = de3211.file.file.Cz(N_start+1:end) - file.Cz_est_val;
file.Cm_residuals = de3211.file.file.Cm(N_start+1:end) - file.Cm_est_val;

% file.Cy_residuals = val_data_asym.file.file.Cy(N_start+1:end) - file.Cy_est_val;
% file.Cl_residuals = val_data_asym.file.file.Cl(N_start+1:end) - file.Cl_est_val;
% file.Cn_residuals = val_data_asym.file.file.Cn(N_start+1:end) - file.Cn_est_val;

file.Cx_mean_res = mean(file.Cx_residuals);
file.Cz_mean_res = mean(file.Cz_residuals);
file.Cm_mean_res = mean(file.Cm_residuals);
% file.Cy_mean_res = mean(file.Cy_residuals);
% file.Cl_mean_res = mean(file.Cl_residuals);
% file.Cn_mean_res = mean(file.Cn_residuals);

save("Training_data/symmetric.mat","file");
fprintf("Done with the symmetric aerdoynamic model estimation\n");
clear("all")
N_start = 200;     % remove first 2 seconds due to bias estimate not being accurate
b = 13.3250;       % [m]
c = 1.9910;        % [m]


%% Asymmetric %%
%%%%%%%%%%%%%%%%
%% Cy, Cl, Cn %%
%%%%%%%%%%%%%%%%

dr3211 = load("Measurement_vector_data/dr3211");
da3211 = load("Measurement_vector_data/da3211");
drdoublet = load("Measurement_vector_data/drdoublet");
dadoublet = load("Measurement_vector_data/dadoublet");

% Concetanate the data of the maneuvers together
% Note that the first 2 seconds fro both maneuvers are discarded due to initial bad bias estimate 
file.Cy_training = [dr3211.file.file.Cy(N_start+1:end);da3211.file.file.Cy(N_start+1:end)];
file.Cl_training = [dr3211.file.file.Cl(N_start+1:end);da3211.file.file.Cl(N_start+1:end)];
file.Cn_training = [dr3211.file.file.Cn(N_start+1:end);da3211.file.file.Cn(N_start+1:end)];

% Add validation data (with drdoublet and dadoublet together) %
file.Cy_val = [drdoublet.file.file.Cy(N_start+1:end);dadoublet.file.file.Cy(N_start+1:end)];
file.Cl_val = [drdoublet.file.file.Cl(N_start+1:end);dadoublet.file.file.Cl(N_start+1:end)];
file.Cn_val = [drdoublet.file.file.Cn(N_start+1:end);dadoublet.file.file.Cn(N_start+1:end)];

beta_training_as = [dr3211.file.file.beta_KF(N_start+1:end-1);da3211.file.file.beta_KF(N_start+1:end-1)];
p_training_as = [dr3211.file.file.pm(N_start+1:end-1);da3211.file.file.pm(N_start+1:end-1)];
V_training_as = [dr3211.file.file.Vm_KF(N_start+1:end-1);da3211.file.file.Vm_KF(N_start+1:end-1)];
r_training_as = [dr3211.file.file.rm(N_start+1:end-1);da3211.file.file.rm(N_start+1:end-1)];
aileron_training_as = [dr3211.file.file.da(N_start+1:end-1);da3211.file.file.da(N_start+1:end-1)];
rudder_training_as = [dr3211.file.file.dr(N_start+1:end-1);da3211.file.file.dr(N_start+1:end-1)];

% VALIDATION
beta_val = [drdoublet.file.file.beta_KF(N_start+1:end-1);dadoublet.file.file.beta_KF(N_start+1:end-1)];
p_val = [drdoublet.file.file.pm(N_start+1:end-1);dadoublet.file.file.pm(N_start+1:end-1)];
V_val = [drdoublet.file.file.Vm_KF(N_start+1:end-1);dadoublet.file.file.Vm_KF(N_start+1:end-1)];
r_val = [drdoublet.file.file.rm(N_start+1:end-1);dadoublet.file.file.rm(N_start+1:end-1)];
aileron_val = [drdoublet.file.file.da(N_start+1:end-1);dadoublet.file.file.da(N_start+1:end-1)];
rudder_val = [drdoublet.file.file.dr(N_start+1:end-1);dadoublet.file.file.dr(N_start+1:end-1)];

% Cy %% 
%theta = [Cy0, CYb, Cyp,Cyr,CYdelta_a,CYdelta_r]
file.A_Cy_training = ones(length(file.Cy_training),6);
file.A_Cy_training(:,2:6) = [beta_training_as,((p_training_as).*b)./(2*V_training_as),(r_training_as.*b)./(2*V_training_as), aileron_training_as,rudder_training_as];

% Cl %%
%theta = [Cl0, Clb, Clp,Clr,Cldelta_a,Cldelta_r]
file.A_Cl_training  = ones(length(file.Cl_training),6);
file.A_Cl_training(:,2:6) = [beta_training_as,((p_training_as).*b)./(2*V_training_as),(r_training_as.*b)./(2*V_training_as), aileron_training_as,rudder_training_as];

% Cn %
%theta = [Cn0, Cnb, Cnp, Cnr, Cndelta_a, Cndelta_r]
file.A_Cn_training  = ones(length(file.Cn_training),6);
file.A_Cn_training(:,2:6) = [beta_training_as,((p_training_as).*b)./(2*V_training_as),(r_training_as.*b)./(2*V_training_as), aileron_training_as,rudder_training_as];

% VALIDATION %
% Cy %% 
%theta = [Cy0, CYb, Cyp,Cyr,CYdelta_a,CYdelta_r]
file.A_Cy_val = ones(length(file.Cy_val),6);
file.A_Cy_val(:,2:6) = [beta_val,((p_val).*b)./(2*V_val),(r_val.*b)./(2*V_val), aileron_val,rudder_val];

% Cl %%
%theta = [Cl0, Clb, Clp,Clr,Cldelta_a,Cldelta_r]
file.A_Cl_val= ones(length(file.Cl_val),6);
file.A_Cl_val(:,2:6) = [beta_val,((p_val).*b)./(2*V_val),(r_val.*b)./(2*V_val), aileron_val,rudder_val];

% Cn %
%theta = [Cn0, Cnb, Cnp, Cnr, Cndelta_a, Cndelta_r]
file.A_Cn_val  = ones(length(file.Cn_val),6);
file.A_Cn_val(:,2:6) = [beta_val,((p_val).*b)./(2*V_val),(r_val.*b)./(2*V_val), aileron_val,rudder_val];


% Parameter estimation %
file.theta_hat_Cy = (file.A_Cy_training' * file.A_Cy_training)^-1 *file.A_Cy_training'*file.Cy_training; 
file.theta_hat_Cl = (file.A_Cl_training' * file.A_Cl_training)^-1 *file.A_Cl_training'*file.Cl_training; 
file.theta_hat_Cn = (file.A_Cn_training' * file.A_Cn_training)^-1 *file.A_Cn_training'*file.Cn_training; 

% training data estimates 
file.Cy_est_training = file.A_Cy_training*file.theta_hat_Cy;
file.Cl_est_training = file.A_Cl_training*file.theta_hat_Cl;
file.Cn_est_training = file.A_Cn_training*file.theta_hat_Cn;

file.Cy_res_training = file.Cy_training - file.Cy_est_training;
file.Cl_res_training = file.Cl_training - file.Cl_est_training;
file.Cn_res_training = file.Cn_training - file.Cn_est_training;

file.Cy_mean_res_training = mean(file.Cy_res_training);
file.Cl_mean_res_training = mean(file.Cl_res_training);
file.Cn_mean_res_training = mean(file.Cn_res_training);

% VALIDATION
file.Cy_est_val = file.A_Cy_val *file.theta_hat_Cy;
file.Cl_est_val  = file.A_Cl_val *file.theta_hat_Cl;
file.Cn_est_val  = file.A_Cn_val *file.theta_hat_Cn;

file.Cy_res_val = file.Cy_val  - file.Cy_est_val ;
file.Cl_res_val  = file.Cl_val  - file.Cl_est_val ;
file.Cn_res_val  = file.Cn_val - file.Cn_est_val ;

file.Cy_mean_res_val  = mean(file.Cy_res_val );
file.Cl_mean_res_val  = mean(file.Cl_res_val );
file.Cn_mean_res_val  = mean(file.Cn_res_val );


save("Training_data/asymmetric.mat","file");
fprintf("Done with the asymmetric aerdoynamic model estimation\n");