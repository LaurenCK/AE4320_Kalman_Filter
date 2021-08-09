%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of the dimensionless aerodynamic force- and moment components 
% Then an OLS regression is performed on the aerodynamic force- and moment
% models as given in the assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Generation of the dimensionless aerodynamic force- and moment components %%
file_names = ["drdoublet.mat", "dr3211.mat","de3211.mat", "da3211.mat","dadoublet.mat"];
Ixx = 11187.8;  % [kg/m^2]
Iyy = 22854.8;
Izz = 31974.8;
Ixz = 1930.1;
b = 13.3250;    % [m]
S = 24.9900;    % [m^2}
c = 1.9910;     % [m]
m = 4500;       % [kg]
z0 = 2000;      % [m] Starting altitude
plotnumber = 0;
for i=1:length(file_names)
    fprintf(strcat("Loading file " ,string(i),"\n"));
    load(strcat("Kalman_data/",file_names(i)));
    altitude = z0 + file.file.z;
    [T, a, P, rho] = atmosisa(altitude);
    rho = rho';
    
    % subtract IMU biases per time step
    file.file.Axm = file.file.Axm - file.file.XX_k1_k1(10,:)';
    file.file.Aym = file.file.Aym - file.file.XX_k1_k1(11,:)';
    file.file.Azm = file.file.Azm - file.file.XX_k1_k1(12,:)';
    file.file.pm = file.file.pm - file.file.XX_k1_k1(13,:)';
    file.file.qm = file.file.qm - file.file.XX_k1_k1(14,:)';
    file.file.rm = file.file.rm - file.file.XX_k1_k1(15,:)';
    
    % subtract IMU biases final value KF
%     file.file.Axm = file.file.Axm - ones(length(file.file.Axm),1).*file.file.XX_k1_k1(10,end)';
%     file.file.Aym = file.file.Aym -ones(length(file.file.Aym),1).*file.file.XX_k1_k1(11,end)';
%     file.file.Azm = file.file.Azm -ones(length(file.file.Azm),1).*file.file.XX_k1_k1(12,end)';
%     file.file.pm = file.file.pm -ones(length(file.file.pm),1).*file.file.XX_k1_k1(13,end)';
%     file.file.qm = file.file.qm -ones(length(file.file.qm),1).*file.file.XX_k1_k1(14,end)';
%     file.file.rm = file.file.rm -ones(length(file.file.rm),1).*file.file.XX_k1_k1(15,end)';
    
    
    % calculate numerical differentiation of p,q,r to get p_dot,q_dot,r_dot
    file.file.p_dot = diff(file.file.pm)/0.01; % !size 12000 instead of 120001
    file.file.q_dot = diff(file.file.qm)/0.01;
    file.file.r_dot = diff(file.file.rm)/0.01;
    
        
    file.file.beta_dot = diff(file.file.betam)/0.01;
    file.file.alpha_dot = diff(file.file.alpham)/0.01;
    file.file.de_dot = diff(file.file.de)/0.01;
    file.file.da_dot = diff(file.file.da)/0.01;
    file.file.dr_dot = diff(file.file.dr)/0.01;
    
    % calculate aerodynamic force- and moment coefficients
    % !!! REMOVE LAST ELEMENT TO MATCH SIZE OF p_dot,r_dot,q_dot !!!
    %h = calc_h(file.file.XX_k1_k1);
    u = file.file.XX_k1_k1(4,:)';
    v = file.file.XX_k1_k1(5,:)';
    w = file.file.XX_k1_k1(6,:)';
    Vm_KF = sqrt(u.*u +v.*v+w.*w);
    file.file.alpha_KF = atan(w./u);
    file.file.beta_KF = atan(v./(sqrt(u.*u+v.*v)));
    file.file.Vm_KF = Vm_KF;
    file.file.Cx = m.*file.file.Axm(1:end-1)./(0.5*rho(1:end-1).*Vm_KF(1:end-1).^2.*S);
    file.file.Cy = m.*file.file.Aym(1:end-1)./(0.5*rho(1:end-1).*Vm_KF(1:end-1).^2.*S);
    file.file.Cz = m.*file.file.Azm(1:end-1)./(0.5*rho(1:end-1).*Vm_KF(1:end-1).^2.*S);
    file.file.u_KF = u;
    file.file.v_KF = v;
    file.file.w_KF = w;
    file.file.u_dot = diff(file.file.u_KF)/0.01;
    file.file.v_dot = diff(file.file.v_KF)/0.01;
    file.file.w_dot = diff(file.file.w_KF)/0.01;  
    
    file.file.Cl = (file.file.p_dot.*Ixx+file.file.qm(1:end-1).*file.file.rm(1:end-1).*(Izz-Iyy)-(file.file.pm(1:end-1).*file.file.qm(1:end-1)+file.file.r_dot).*Ixz)./(0.5*rho(1:end-1).*Vm_KF(1:end-1).^2.*S.*b);
    file.file.Cm = (file.file.q_dot.*Iyy+file.file.rm(1:end-1).*file.file.pm(1:end-1).*(Ixx-Izz)+(file.file.pm(1:end-1).^2-file.file.rm(1:end-1).^2).*Ixz)./(0.5*rho(1:end-1).*Vm_KF(1:end-1).^2.*S.*c);
    file.file.Cn = (file.file.r_dot.*Izz+file.file.pm(1:end-1).*file.file.qm(1:end-1).*(Iyy-Ixx)+(file.file.qm(1:end-1).*file.file.rm(1:end-1)-file.file.p_dot).*Ixz)./(0.5*rho(1:end-1).*Vm_KF(1:end-1).^2.*S.*b);
    
    fprintf(strcat("Saving file " ,string(i),"\n"));
    save(strcat("Measurement_vector_data\",file_names{i}),"file");
    
%     plotnumber = plotnumber + 1;
%     figure(plotnumber)
%     plot(file.file.t(1:end-1),file.file.Cx);
%     hold on 
%     plot(file.file.t(1:end-1),file.file.Cz);
%     hold on
%     plot(file.file.t(1:end-1),file.file.Cm);
%     hold on
%     plot(file.file.t(1:end-1),file.file.Cy);
%     hold on
%     plot(file.file.t(1:end-1),file.file.Cl);
%     hold on
%     plot(file.file.t(1:end-1),file.file.Cn);
%     title(strcat("Input data: ",file_names(i)));
%     legend("Cx","Cz","Cm","Cy","Cl","Cn");
    
end
