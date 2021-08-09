%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Calculate parameter covariance matrices, print and plot %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear("all")
close("all")

asym = load("Training_data/asymmetric.mat"); 
sym =  load("Training_data/symmetric.mat"); 


%% Parameter covariance matrix %%
% covariance of the parameters are on the diagonal %
% cross-covariance between parameters are the off-diagonal terms %

P_Cx = (sym.file.A_Cx_val'*sym.file.A_Cx_val)^-1;
P_Cz = (sym.file.A_Cz_val'*sym.file.A_Cz_val)^-1;
P_Cm = (sym.file.A_Cm_val'*sym.file.A_Cm_val)^-1;

P_Cy = (asym.file.A_Cy_val'*asym.file.A_Cy_val)^-1;
P_Cl = (asym.file.A_Cl_val'*asym.file.A_Cl_val)^-1;
P_Cn = (asym.file.A_Cn_val'*asym.file.A_Cn_val)^-1;

var_Cx_res = var(sym.file.Cx_residuals);
var_Cz_res = var(sym.file.Cz_residuals);
var_Cm_res = var(sym.file.Cm_residuals);

var_Cy_res = var(asym.file.Cy_res_val);
var_Cl_res = var(asym.file.Cl_res_val);
var_Cn_res = var(asym.file.Cn_res_val);

Cx_cov = var_Cx_res*P_Cx
Cz_cov = var_Cz_res*P_Cz
Cm_cov = var_Cm_res*P_Cm

Cy_cov = var_Cy_res*P_Cy
Cl_cov = var_Cl_res*P_Cl
Cn_cov = var_Cn_res*P_Cn


figure(1)
sgtitle("Parameter Covariance Matrices")
subplot(2,3,1)
imagesc(Cx_cov)
colorbar
title("Cx")
hold on
subplot(2,3,2)
imagesc(Cz_cov)
colorbar
title("Cz")
hold on
subplot(2,3,3)
imagesc(Cm_cov)
colorbar
title("Cm")
hold on
subplot(2,3,4)
imagesc(Cy_cov)
colorbar
title("Cy")
hold on
subplot(2,3,5)
imagesc(Cl_cov)
colorbar
title("Cl")
hold on
subplot(2,3,6)
imagesc(Cn_cov)
colorbar
title("Cn")



