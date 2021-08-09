
function data_generation(sensitivity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of data to be used in the state estimation (Kalman Filter).
% Noise and bias terms are added. Data is saved in Data and Data2 ( = 10x
% higher airdata noise level when sensitivity == 1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_names = ["drdoublet.mat", "dr3211.mat","de3211.mat", "da3211.mat","dadoublet.mat"];
rng('default'); % seed
plotnumber = 10;

% wind speed [m/s]
Wx =  -10; 
Wy = 3; 
Wz = 2;

% Real biases of accelerometer and rate gyros
gamma_xr = 0.01;                % [m/s^2]
gamma_yr = 0.01;
gamma_zr = 0.01;
gamma_pr = 0.001*pi/180;        % [rad/s] !
gamma_qr = 0.001*pi/180;  
gamma_rr = 0.001*pi/180; 

% Noise standard deviation of accelerometer and rate gyros
sigma_xr = 0.01;                % [m/s^2]
sigma_yr = 0.01;
sigma_zr = 0.01;
sigma_pr = 0.001*pi/180;        % [rad/s] !
sigma_qr = 0.001*pi/180;  
sigma_rr = 0.001*pi/180; 

% GPS noise estandard deviations
sigma_GPS_pos = 5;              % [m]
sigma_GPS_vel = 0.1;            % [m/s]
sigma_GPS_att = 0.1*pi/180;     % [rad] !

% Airdata standard deviations
if sensitivity==1
    sigma_airdata_vtas = 10*0.1;       % Vtas [m/s]
    sigma_airdata_angle = 10*0.1*pi/180; %AoA and beta [rad], [rad]!
else
    sigma_airdata_vtas = 0.1;       % Vtas [m/s]
    sigma_airdata_angle = 0.1*pi/180; %AoA and beta [rad], [rad]!
end
    
for i = 1:length(file_names)
   file = load(file_names(i));
   
   % Integrating groundspeed in inertial frame to aircraft position
   % (includes wind)
   file.x = cumtrapz(file.t,file.u_n+Wx) ;
   file.y = cumtrapz(file.t,file.v_n+Wy) ;
   file.z = cumtrapz(file.t,file.w_n+Wz) ;
   
   % Generating IMU data including a BIAS and NOISE
   file.Axm = file.Ax + gamma_xr + normrnd(0,sigma_xr, size(file.Ax));
   file.Aym = file.Ay + gamma_yr + normrnd(0,sigma_yr,size(file.Ay));
   file.Azm = file.Az + gamma_zr + normrnd(0,sigma_zr,size(file.Az));
   
   file.pm = file.p + gamma_pr + normrnd(0,sigma_pr,size(file.p));
   file.qm = file.q + gamma_qr + normrnd(0,sigma_qr,size(file.q));
   file.rm = file.r + gamma_rr + normrnd(0,sigma_rr,size(file.r));
   
   % Generating GPS data including NOISE
   file.xGPSm = file.x + normrnd(0,sigma_GPS_pos,size(file.x));
   file.yGPSm = file.y + normrnd(0,sigma_GPS_pos,size(file.y));
   file.zGPSm = file.z + normrnd(0,sigma_GPS_pos,size(file.z));
   
   file.uGPSm = file.u_n + Wx + normrnd(0,sigma_GPS_vel,size(file.u_n));
   file.vGPSm = file.v_n + Wy + normrnd(0,sigma_GPS_vel,size(file.v_n));
   file.wGPSm = file.w_n + Wz + normrnd(0,sigma_GPS_vel,size(file.w_n));
   
   file.phiGPSm = file.phi + normrnd(0,sigma_GPS_att,size(file.phi));
   file.thetaGPSm = file.theta + normrnd(0,sigma_GPS_att,size(file.theta));
   file.psiGPSm = file.psi + normrnd(0,sigma_GPS_att,size(file.psi));
   
   % Generating the airdata including NOISE
   file.Vm = file.vtas + normrnd(0,sigma_airdata_vtas,size(file.vtas));
   file.alpham = file.alpha + normrnd(0,sigma_airdata_angle,size(file.alpha));
   file.betam = file.beta + normrnd(0,sigma_airdata_angle,size(file.beta));
   
   file.real_states = [file.x';file.y';file.z';file.u_n';file.v_n';file.w_n';file.phi';file.theta';file.psi'];
    % u,v,w inertial
   if sensitivity == 0
        save(strcat("Data\",file_names{i}),"file");
   end
   
   if sensitivity == 1
       save(strcat("Data2\",file_names{i}),"file");
   end
   
   
   if i == 4
    close("all")
    plotnumber = plotnumber +1;
    figure(plotnumber)
    sgtitle("Aircraft Position Generated from Simulated Airspeed Data and Wind Velocity for da3211.mat")
    subplot(1,2,1)
    plot(file.t,file.u_n+Wx,"--")
    hold on
    plot(file.t,file.v_n+Wy,"--")
    hold on
    plot(file.t,file.w_n+Wz,"--")
    hold on
    plot(file.t, file.x)
    hold on
    plot(file.t, file.y)
    hold on
    plot(file.t, file.z)
    %legend("u_{ground}","v_{ground}","w_{ground}","x","y","z");     
    xlabel("Time [s]")
    grid on
    title("Position and Ground Velocities")
    
    subplot(1,2,2)    
    plot(file.t,file.u_n+Wx,"--")
    hold on
    plot(file.t,file.v_n+Wy,"--")
    hold on
    plot(file.t,file.w_n+Wz,"--")
    hold on
    plot(file.t, file.x)
    hold on
    plot(file.t, file.y)
    hold on
    plot(file.t, file.z)
    title("Zoomed-in on Ground Velocities")
    legend("u_{ground}","v_{ground}","w_{ground}","x","y","z");     
    xlabel("Time [s]")
    ylim([-100 100])
    grid on
    
   end
   
end


    
%end function
end



















    



