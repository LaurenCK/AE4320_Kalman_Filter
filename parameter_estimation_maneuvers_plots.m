%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot all results from the parameter estimation trained on the training data
% Using both training and the validation data. 
% Calculate residuals and autocorrelation functions.
% Code is seperated into  1) symmetric and 2) asymmetric part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close("all");
clear("all")
N_start = 200;     % remove first 2 seconds due to bias estimate not being accurate
b = 13.3250;       % [m]
c = 1.9910;        % [m]
plotnumber = 0;

%%
asym = load("Training_data/asymmetric.mat"); % contains theta_hats, A_train and A_valid, C_estimate: training and validation data
sym =  load("Training_data/symmetric.mat"); % also contains MEASURED C_training and C_valid


%% Plotting TRAINING data estimate CY, Cl, Cn: dr3211 and da3211 file together %%
% Cutting the data to make the plot more clear %
Cy_training = [asym.file.Cy_training(5/0.01:25/0.01);asym.file.Cy_training(130/0.01:end-9500)];
Cl_training = [asym.file.Cl_training(5/0.01:25/0.01);asym.file.Cl_training(130/0.01:end-9500)];
Cn_training = [asym.file.Cn_training(5/0.01:25/0.01);asym.file.Cn_training(130/0.01:end-9500)];

Cy_est_training = [asym.file.Cy_est_training(5/0.01:25/0.01);asym.file.Cy_est_training(130/0.01:end-9500)];
Cl_est_training = [asym.file.Cl_est_training(5/0.01:25/0.01);asym.file.Cl_est_training(130/0.01:end-9500)];
Cn_est_training = [asym.file.Cn_est_training(5/0.01:25/0.01);asym.file.Cn_est_training(130/0.01:end-9500)];

t_trimmed = (0:0.01:(length(Cy_training)-1)*0.01)';

plotnumber = plotnumber + 1;
figure(plotnumber)
sgtitle("Estimates on Training Data dr3211.mat and da3211.mat");
subplot(3,1,1)
plot(t_trimmed,Cy_training,"--")
hold on
plot(t_trimmed,Cy_est_training)
title("Cy")
legend("Cy","Cy_{estimate}");
grid on
xlabel("Time [s]")
xlim([0 32])
ylabel("[-]")
%xlim([8 30])
%ylim([-0.03,0.03])
hold on
subplot(3,1,2)
plot(t_trimmed,Cl_training,"--")
hold on
plot(t_trimmed,Cl_est_training)
legend("Cl","Cl_{estimate}");
grid on
title("Cl")
xlabel("Time [s]")
xlim([0, 32])
ylabel("[-]")
%xlim([8 30])
hold on
subplot(3,1,3)
plot(t_trimmed,Cn_training,"--")
hold on
plot(t_trimmed,Cn_est_training)
legend("Cn","Cn_{estimate}");
grid on
title("Cn")
xlabel("Time [s]")
xlim([0, 32])
ylabel("[-]")
%xlim([8 30])


%% ----- plotting for SYMMETRIC (training&validating data) ------ %
data = load("Data/drdoublet.mat");
t_sym = data.file.t(N_start+1:end-1);
plotnumber = plotnumber + 1;
figure(plotnumber)
sgtitle("Estimates on Training/Validation Data de3211.mat");
subplot(3,1,1)
plot(t_sym ,sym.file.Cx_val,"--")
hold on
plot(t_sym ,sym.file.Cx_est_val)
title("Cx")
legend("Cx","Cx_{estimate}");
grid on
xlabel("Time [s]")
ylabel("[-]")
%xlim([8 30])
%ylim([-0.03,0.03])
hold on
subplot(3,1,2)
plot(t_sym ,sym.file.Cz_val,"--")
hold on
plot(t_sym ,sym.file.Cz_est_val)
legend("Cz","Cz_{estimate}");
grid on
title("Cz")
xlabel("Time [s]")
ylabel("[-]")
%xlim([8 30])
hold on
subplot(3,1,3)
plot(t_sym ,sym.file.Cm_val,"--")
hold on
plot(t_sym ,sym.file.Cm_est_val)
legend("Cm","Cm_{estimate}");
grid on
title("Cm")
xlabel("Time [s]")
ylabel("[-]")
%xlim([8 30])


%% plotting RESIDUAL for TOTAL training data dr3211 and da3211. AND de3211.mat %%

plotnumber = plotnumber + 1; 
figure(plotnumber)
sgtitle("Model Residuals on Training Data")
subplot(2,3,1)
plot(t_sym,sym.file.Cx_residuals)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(2),ylim(2),strcat("mean: ",string(sym.file.Cx_mean_res)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cx")
xlabel("Time [s]")
ylabel("[-]")
hold on
subplot(2,3,2)
plot(t_sym,sym.file.Cz_residuals)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(2),ylim(2),strcat("mean: ",string(sym.file.Cz_mean_res)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cz")
xlabel("Time [s]")
ylabel("[-]")
hold on
subplot(2,3,3)
plot(t_sym,sym.file.Cm_residuals)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(2),ylim(2),strcat("mean: ",string(sym.file.Cm_mean_res)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cm")
xlabel("Time [s]")
ylabel("[-]")
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_asym = (0:0.01:(length(asym.file.Cy_training)-1)*0.01)';
subplot(2,3,4)
plot(t_asym,asym.file.Cy_res_training)
ylim_axis=get(gca,'ylim');
xlim_axis=get(gca,'xlim');
text(xlim_axis(2)*1.06,ylim_axis(2),strcat("mean: ",string(asym.file.Cy_mean_res_training)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cy")
xlabel("Time [s]")
ylabel("[-]")
hold on
subplot(2,3,5)
plot(t_asym,asym.file.Cl_res_training)
ylim_axis=get(gca,'ylim');
xlim_axis=get(gca,'xlim');
text(xlim_axis(2)*1.06,ylim_axis(2)*1.3,strcat("mean: ",string(asym.file.Cl_mean_res_training)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cl")
xlabel("Time [s]")
ylabel("[-]")
hold on
subplot(2,3,6)
plot(t_asym,asym.file.Cn_res_training)
ylim_axis=get(gca,'ylim');
xlim_axis=get(gca,'xlim');
text(xlim_axis(2)*1.06,ylim_axis(2),strcat("mean: ",string(asym.file.Cn_mean_res_training)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cn")
xlabel("Time [s]")
ylabel("[-]")
clear("xlim","ylim");


%% Check autocorrelation of TRAINING residuals %%
% xcorr
% 
conf = ones(2*11800+1,1).*1.96/sqrt(length(sym.file.Cx_residuals));
plotnumber = plotnumber+1;
    figure(plotnumber)
    sgtitle("Cross-correlation of the Model Residuals on Training Data")
    subplot(2,3,1)
    [lags,c] = xcorr(sym.file.Cx_residuals,11800,"normalized");
    plot(c,lags);
    hold on
    plot(c,conf,"r--")
    hold on
    plot(c,-conf,"r--")
    title("Cx")
    grid on
    xlabel("Number of Lags")
    ylabel("[-]")
    hold on
    subplot(2,3,2)
    [lags,c] = xcorr(sym.file.Cz_residuals,11800,"normalized");
    plot(c,lags);
    hold on
    plot(c,conf,"r--")
    hold on
    plot(c,-conf,"r--")
    title("Cz")
    grid on
    xlabel("Number of Lags")
    ylabel("[-]")
    hold on
    subplot(2,3,3)
    [lags,c] = xcorr(sym.file.Cm_residuals,11800,"normalized");
    plot(c,lags);
    hold on
    plot(c,conf,"r--")
    hold on
    plot(c,-conf,"r--")
    title("Cm")
    grid on
    xlabel("Number of Lags")
    ylabel("[-]")
    hold on
    subplot(2,3,4)
    [lags,c] = xcorr(asym.file.Cy_res_training,11800,"normalized");
    plot(c,lags);
    hold on
    plot(c,conf,"r--")
    hold on
    plot(c,-conf,"r--")
    title("Cy")
    grid on
    xlabel("Number of Lags")
    ylabel("[-]")
    hold on
    subplot(2,3,5)
    [lags,c] = xcorr(asym.file.Cl_res_training,11800,"normalized");
    plot(c,lags);
    hold on
    plot(c,conf,"r--")
    hold on
    plot(c,-conf,"r--")
    title("Cl")
    grid on
    xlabel("Number of Lags")
    ylabel("[-]")
    hold on 
    subplot(2,3,6)
    [lags,c] = xcorr(asym.file.Cn_res_training,11800,"normalized");
    plot(c,lags);
    hold on
    plot(c,conf,"r--")
    hold on
    plot(c,-conf,"r--")
    title("Cn")
    grid on
    xlabel("Number of Lags")
    ylabel("[-]")
    legend("Cross-correlation","95% confidence bounds")
    

%%%%%%%%%%%%%%%
%% VALIDATION%%
%%%%%%%%%%%%%%%

%% Plotting VALIDATION data estimate CY, Cl, Cn: drdoublet and dadoublet file together %%
% Cutting the data to make the plot more clear %
Cy_val = [asym.file.Cy_val(5/0.01:20/0.01);asym.file.Cy_val(125/0.01:end-10000)];
Cl_val = [asym.file.Cl_val(5/0.01:20/0.01);asym.file.Cl_val(125/0.01:end-10000)];
Cn_val = [asym.file.Cn_val(5/0.01:20/0.01);asym.file.Cn_val(125/0.01:end-10000)];

Cy_est_val = [asym.file.Cy_est_val(5/0.01:20/0.01);asym.file.Cy_est_val(125/0.01:end-10000)];
Cl_est_val = [asym.file.Cl_est_val(5/0.01:20/0.01);asym.file.Cl_est_val(125/0.01:end-10000)];
Cn_est_val = [asym.file.Cn_est_val(5/0.01:20/0.01);asym.file.Cn_est_val(125/0.01:end-10000)];

t_trimmed = (0:0.01:(length(Cy_val)-1)*0.01)';

plotnumber = plotnumber + 1;
figure(plotnumber)
sgtitle("Estimates on Validation Data drdoublet.mat and dadoublet.mat");
subplot(3,1,1)
plot(t_trimmed,Cy_val,"--")
hold on
plot(t_trimmed,Cy_est_val)
title("Cy")
legend("Cy","Cy_{estimate}");
grid on
xlabel("Time [s]")
xlim([0 27])
ylabel("[-]")
%xlim([8 30])
%ylim([-0.03,0.03])
hold on
subplot(3,1,2)
plot(t_trimmed,Cl_val,"--")
hold on
plot(t_trimmed,Cl_est_val)
legend("Cl","Cl_{estimate}");
grid on
title("Cl")
xlabel("Time [s]")
xlim([0 27])
ylabel("[-]")
%xlim([8 30])
hold on
subplot(3,1,3)
plot(t_trimmed,Cn_val,"--")
hold on
plot(t_trimmed,Cn_est_val)
legend("Cn","Cn_{estimate}");
grid on
title("Cn")
xlabel("Time [s]")
xlim([0 27])
ylabel("[-]")
%xlim([8 30])

%% plotting RESIDUAL for TOTAL validation data dr3211 and da3211. NO de3211.mat %%

plotnumber = plotnumber + 1; 
figure(plotnumber)
sgtitle("Model Residuals on Validation Data drdoublet.mat and dadoubletmat");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_asym = (0:0.01:(length(asym.file.Cy_training)-1)*0.01)';
subplot(3,1,1)
plot(t_asym,asym.file.Cy_res_val)
ylim_axis=get(gca,'ylim');
xlim_axis=get(gca,'xlim');
text(xlim_axis(2),ylim_axis(2),strcat("mean: ",string(asym.file.Cy_mean_res_val)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cy")
xlabel("Time [s]")
ylabel("[-]")
hold on
subplot(3,1,2)
plot(t_asym,asym.file.Cl_res_val)
ylim_axis=get(gca,'ylim');
xlim_axis=get(gca,'xlim');
text(xlim_axis(2),ylim_axis(2)*0.95,strcat("mean: ",string(asym.file.Cl_mean_res_val)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cl")
xlabel("Time [s]")
ylabel("[-]")
hold on
subplot(3,1,3)
plot(t_asym,asym.file.Cn_res_val)
ylim_axis=get(gca,'ylim');
xlim_axis=get(gca,'xlim');
text(xlim_axis(2),ylim_axis(2),strcat("mean: ",string(asym.file.Cn_mean_res_val)),'VerticalAlignment','top', 'HorizontalAlignment','right')
grid on
title("Cn")
xlabel("Time [s]")
ylabel("[-]")
clear("xlim","ylim");

theta_hat_final = [string(sym.file.theta_hat_Cx);string(sym.file.theta_hat_Cz);string(sym.file.theta_hat_Cm);string(asym.file.theta_hat_Cy);string(asym.file.theta_hat_Cl);string(asym.file.theta_hat_Cn)];
parameter_names = ["Cx0","Cxalpha","Cxalpha^2","Cxq","Cxdelta_e","CxTc","Cz0","Czalpha","Czq","Czdelta_e","CzTc","Cm0","Cmalpha","Cmq","Cmdelta_e","CmTc","Cy0","Cybeta","Cyp","Cyr","Cydelta_a","Cy_delta_r","Cl0","Clbeta","Clp","Clr","Cldelta_a","Cl_delta_r","Cn0","Cnbeta","Cnp","Cnr","Cndelta_a","Cn_delta_r"];
table_theta_hats = table(parameter_names',theta_hat_final);
display(table_theta_hats)
%table2latex(table_theta_hats, 'Training_data/table_hats.txt')