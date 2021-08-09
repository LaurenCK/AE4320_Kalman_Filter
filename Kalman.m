

function Kalman(std_VTAS,std_alpha, std_beta,sensitivity)
%% Extended- and iterated Kalman Filter script %%%%%%%%%%
% When sensitivity == 1 the data is saved into Data2%%%
% Most parts taken from AE4320 course
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n               = 18;           % number of states
nm              = 12;           % number of measurements
%m               = 6;            % number of inputs
%dt              = 0.01;        % time step (s)
%N               = 1000;        % sample size
epsilon         = 1e-10;
runIEKF          = 1;           % set 1 for IEKF and 0 for EKF
maxIterations   = 10;

file_names = ["drdoublet.mat", "dr3211.mat","de3211.mat", "da3211.mat","dadoublet.mat"];

%% Iterate over each maneuver and for each maneuver do a state estimation on all 18 states %%
for i = 1:length(file_names)
    
   if sensitivity == 0
       try
           file = load(strcat("Data/",file_names{i}));
       catch
           fprintf("No data files found");
       end
   end
   
   if sensitivity == 1
       try
           file = load(strcat("Data2/",file_names{i}));
       catch
           fprintf("No data files found");
       end
   end
   
           
   N = length(file.file.t);    % number of data points
   dt = file.file.t(2)-file.file.t(1);
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set initial values for states and statistics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E_x_0       = [file.file.xGPSm(1); % initial state for all 18 states, use initial values of GPS
                   file.file.yGPSm(1);   
                   file.file.zGPSm(1);  
                   file.file.uGPSm(1);   
                   file.file.vGPSm(1);   
                   file.file.wGPSm(1);   
                   0;%file.file.phiGPSm(1);
                   0;%file.file.thetaGPSm(1);
                   0;%file.file.psiGPSm(1);
                   0;0;0;0;0;0;0;0;0];
    
    %%E_x_0 = 0*ones(18,1);

    %G           = calc_G(E_x_0);  % noise input matrix at x_0

    % Initial estimate for estimation error covariance matrix
    std_x_0     = 1*ones(1,18);
    std_x_0(10:18) = 10;%10^-6;
    P_0         = diag(std_x_0.^2,0); 

    % System noise statistics (what we guess it is):
    std_x_IMU       = 0.01;  % [m/s*s]
    std_y_IMU       = 0.01;    
    std_z_IMU       = 0.01;    
    std_p_IMU       = 0.001*pi/180; %[rad/s]!  
    std_q_IMU       = 0.001*pi/180;    
    std_r_IMU       = 0.001*pi/180; 
    std_system = [std_x_IMU,std_y_IMU,std_z_IMU,std_p_IMU,std_q_IMU,std_r_IMU];

    Q = diag(std_system.^2,0);      % System noise covariance matrix
 

    % Measurement/sensor noise statistics (what we guess it is):
    std_GPS_x = 5;                  %[m]
    std_GPS_y = 5; 
    std_GPS_z = 5; 
    std_GPS_u = 0.1;                %[m/s]
    std_GPS_v = 0.1; 
    std_GPS_w = 0.1; 
    std_GPS_phi = 0.1*pi/180;       %[rad]!
    std_GPS_theta = 0.1*pi/180;    
    std_GPS_psi = 0.1*pi/180; 
    %std_VTAS = 0.1;                 %[m/s]
    %std_alpha = 0.1*pi/180;         %[rad]!
    %std_beta = 0.1*pi/180;
    std_sensor = [std_GPS_x,std_GPS_y,std_GPS_z,std_GPS_u,std_GPS_v,std_GPS_w, std_GPS_phi,std_GPS_theta,std_GPS_psi,std_VTAS,std_alpha,std_beta];
    
    R = diag(std_sensor.^2,0);      % sensor noise covariance matrix
    

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Iterated Extended Kalman filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t_k             = 0; 
    t_k1            = dt;

    XX_k1_k1        = zeros(n, N);      % Our state estimates
    PP_k1_k1        = zeros(n, n, N);   % Our optimal prediction error covariance
    STD_x_cor       = zeros(n, N);      % Standard deviation of the optimal prediction error covariance matrix
    ZZ_pred         = zeros(nm, N);     % prediction of observation
    STD_z           = zeros(nm,N);
    IEKFitcount     = zeros(N, 1);

    x_k1_k1         = E_x_0;    % x(0|0)=E{x_0}
    P_k1_k1         = P_0;      % P(0|0)=P(0)
    
    % state input at all time steps for the IMU measured Axm, Aym, Azm, pm,qm,rm
    U_k = [file.file.Axm';
           file.file.Aym';
           file.file.Azm';
           file.file.pm';
           file.file.qm';
           file.file.rm'];
    % Obervation equation at all time steps GPS, VTAS, alpha beta
    Z_k = [file.file.xGPSm';
           file.file.yGPSm';
           file.file.zGPSm';
           file.file.uGPSm';
           file.file.vGPSm';
           file.file.wGPSm';
           file.file.phiGPSm';
           file.file.thetaGPSm';
           file.file.psiGPSm';
           file.file.Vm';
           file.file.alpham';
           file.file.betam'];           


%% Run the filter through all N samples %%
for k = 1:N
    % 1) x(k+1|k) (prediction) by integrating the nonlinear state transition function f(x,u,t)
    [t, x_k1_k]     = rk4(@calc_f, x_k1_k1,U_k(:,k), [t_k t_k1]); 
                        %ode45(@calc_f,[t_k t_k1],x_k1_k1); % integrate from t_k to t_k1 at x_k1_k1 (previous optimal state estimation)
    
    % z(k+1|k) (predicted observation)
    z_k1_k          = calc_h(x_k1_k)';           % prediction of observation at the integrated predicted state
    ZZ_pred(:,k)    = z_k1_k;                    % store this observation prediction, since later prediction observations
                                                 % are "corrected" using the actual observation
        
    % 2) Calc Jacobians Phi(k+1,k) and Gamma(k+1, k)
    Fx              = calc_dF(x_k1_k, U_k(:,k)); % perturbation of f(x,u,t)
    
    % 3) discritize PERTUBATION transition- and noise input matrix
    G = calc_G(x_k1_k);
    [Phi, Gamma] = c2d(Fx,G,dt);
    
    % 4) P(k+1|k) (prediction error covariance matrix)
    P_k1_k          = Phi*P_k1_k1*Phi' + Gamma*Q*Gamma'; 
    
    %% Run IEKF %%
    % IEKF does some extra iterations to improve lineatization around an
    % improved nominal state of the observation equation
        
    % Goal: update the nominal states using the current MEASUREMENT and
    % RE-LINEARIZING the system using the IMPROVED NOMINAL STATES.
    if (runIEKF==1)
        
        eta2 = x_k1_k;
        err = 2*epsilon;        
        iterations = 0;
        
        while(err > epsilon)
            if (iterations > maxIterations)
                fprintf('Terminating IEKF: exceeded max iterations (%d)\n', maxIterations);
                break
            end
            iterations = iterations + 1;
            eta1 = eta2; 
                        
            % Construct the Jacobian H = d/dx(h(x))) with h(x) the observation model transition matrix
            Hx      = calc_dH(eta1);   

            % Check observability of state
            if (k == 1 && iterations == 1)
                rankHF = kf_calcObsRank(Hx, Fx);
                if (rankHF < n)
                    warning('The current state is not observable; rank of Observability Matrix is %d, should be %d', rankHF, n);
                end
            end
                    
            
            % P_zz(k+1|k) (covariance matrix of innovation)
            P_zz        = (Hx*P_k1_k * Hx' + R);            % covariance matrix of observation error
            std_z       = sqrt(diag(P_zz));                 % standard deviation of observation error (for validation) 

            % calculate the Kalman gain matrix
            K       = P_k1_k*Hx'/(P_zz);
            % new observation state
            z_p     = calc_h(eta1)';

            eta2    = x_k1_k + K*(Z_k(:,k) - z_p - Hx*(x_k1_k - eta1));
            err     = norm((eta2 - eta1), inf) / norm(eta1, inf);
            
            
        end
      IEKFitcount(k)  = iterations;
      x_k1_k1         = eta2;  
        
    end
    
    %% Run EKF %%
    if (runIEKF==0)
        Hx              = calc_dH(x_k1_k);  % perturbation of h(x,u,t)
        
        % 5) Kalman gain calculation

        % P_zz(k+1|k) (covariance matrix of innovation)
        P_zz            = (Hx*P_k1_k * Hx' + R);           % covariance matrix of observation error
        std_z           = sqrt(diag(P_zz));                 % standard deviation of observation error (for validation)    

        % K(k+1) (gain)
        K               = (P_k1_k * Hx')/(P_zz);
        
        % 6) Measurement update
        % Calculate optimal state x(k+1|k+1) 
        x_k1_k1         = x_k1_k + K*(Z_k(:,k) - z_k1_k);     

        
    end
    
    % P(k+1|k+1) (correction) using the numerically stable form of P_k_1k_1 = (eye(n) - K*Hx) * P_kk_1;
    P_k1_k1         = (eye(n) - K*Hx) * P_k1_k*(eye(n) - K*Hx)' + K*R*K';  
    std_x_cor       = sqrt(diag(P_k1_k1));

    % Next step
    t_k             = t_k1; 
    t_k1            = t_k1 + dt;
    
    % store results
    XX_k1_k1(:,k)     = x_k1_k1;
    PP_k1_k1(:,:,k)   = P_k1_k1;
    STD_x_cor(:,k)    = std_x_cor; 
    STD_z(:,k)        = std_z;
   
   
end
fprintf(strcat("Done with maneuver ",string(i)," saving states...\n"));

%% Save the Estimated states ... for the maneuver %% 

file.file.XX_k1_k1 = XX_k1_k1;
file.file.PP_k1_k1 = PP_k1_k1;
file.file.STD_x_cor = STD_x_cor;


% Measurement estimation error (Innovation) = difference between measured states and KF states
file.file.innovation = ZZ_pred - Z_k;  


% converting the u,v,w body velocities to intertial frame
file.file.u_inertial_KF = (XX_k1_k1(4,:).*cos(XX_k1_k1(8,:))+(XX_k1_k1(5,:).*sin(XX_k1_k1(7,:))+XX_k1_k1(6,:).*cos(XX_k1_k1(7,:))).*sin(XX_k1_k1(8,:))).*cos(XX_k1_k1(9,:))-(XX_k1_k1(5,:).*cos(XX_k1_k1(7,:))-XX_k1_k1(6,:).*sin(XX_k1_k1(7,:))).*sin(XX_k1_k1(9,:));%+XX_k1_k1(16,:);
file.file.v_inertial_KF = (XX_k1_k1(4,:).*cos(XX_k1_k1(8,:))+(XX_k1_k1(5,:).*sin(XX_k1_k1(7,:))+XX_k1_k1(6,:).*cos(XX_k1_k1(7,:))).*sin(XX_k1_k1(8,:))).*sin(XX_k1_k1(9,:))+(XX_k1_k1(5,:).*cos(XX_k1_k1(7,:))-XX_k1_k1(6,:).*sin(XX_k1_k1(7,:))).*cos(XX_k1_k1(9,:));%+XX_k1_k1(17,:);
file.file.w_inertial_KF = -XX_k1_k1(4,:).*sin(XX_k1_k1(8,:))+(XX_k1_k1(5,:).*sin(XX_k1_k1(7,:))+XX_k1_k1(6,:).*cos(XX_k1_k1(7,:))).*cos(XX_k1_k1(8,:));%+XX_k1_k1(18,:);

% State estimation error (REAL error!) = Difference between the KF states and the real states (no bias, no noise) 
file.file.real_error = XX_k1_k1(1:9,:) - file.file.real_states;
file.file.real_error(4,:) = file.file.u_inertial_KF - file.file.real_states(4,:); %u inertial 
file.file.real_error(5,:) = file.file.v_inertial_KF - file.file.real_states(5,:); %v 
file.file.real_error(6,:) = file.file.w_inertial_KF - file.file.real_states(6,:); %w
%file.file.real_error_norm = 2*file.file.real_error./(XX_k1_k1(1:9,:)+file.file.real_states);
file.file.real_error_norm = file.file.real_error./(file.file.real_states+0.0001);


if sensitivity == 1
    save(strcat("Kalman_data2\",file_names{i}),"file");
else
    save(strcat("Kalman_data\",file_names{i}),"file");
end


fprintf("Perfoming KF on next maneuver...\n");


end

%end function
end


