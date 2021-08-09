%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot all KF result plots %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close("all")
file_names = ["drdoublet.mat", "dr3211.mat","de3211.mat", "da3211.mat","dadoublet.mat"];
state_names = ["x","y","z","u","v","w","Phi","Theta","Psi"];
plotnumber = 0;
%%
%data_generation(0);
%Kalman(0.1,0.1*pi/180,0.1*pi/180,0);
%%
for i = 1:1%length(file_names)
     try
        load(strcat("Kalman_data\",file_names{i}));
     catch
         fprintf("File not found, running Kalman filter...");
         run(Kalman_filter.m);
     end
    
%% Kalman Filter results: KF estimate, real, and observation (GPS data)
    plotnumber = plotnumber + 1;
    figure(plotnumber)
    sgtitle(strcat("Kalman Filter Results of ",file_names{i}));

    subplot(3,3,1);
    plot(file.file.t, file.file.xGPSm, "LineWidth",0.5,"color",[0.9290, 0.6940, 0.1250]	); %% intertial!!! ground speed
    hold on   
    plot(file.file.t, file.file.XX_k1_k1(1,:),"color","b");
    hold on
    plot(file.file.t,file.file.x, "--","color","r"); %%  intertial!!! ground speed
    title("Position x");
    xlabel("Time [s]");
    ylabel("Distance [m]");
    xlim([0 120])
    %legend("Real x","Estimated x", "GPS x");

    subplot(3,3,2);
    plot(file.file.t, file.file.yGPSm, "LineWidth",0.5,"color",[0.9290, 0.6940, 0.1250]	);
    hold on
    plot(file.file.t, file.file.XX_k1_k1(2,:),"color","b");  
    hold on
    plot(file.file.t,file.file.y, "--","color","r");
    title("Position y");
    xlabel("Time [s]");
    ylabel("Distance [m]");
    xlim([0 120])
    %legend("Real y","Estimated y", "GPS y");
    
    subplot(3,3,3);
    plot(file.file.t, file.file.zGPSm, "LineWidth",0.5,"color",[0.9290, 0.6940, 0.1250]	);
    hold on
    plot(file.file.t, file.file.XX_k1_k1(3,:),"color","b");
    hold on
    plot(file.file.t,file.file.z, "--","color","r");
    title("Position z");
    xlabel("Time [s]");
    ylabel("Distance [m]");
    xlim([0 120])
    %legend("Real z","Estimated z", "GPS z");
     
    subplot(3,3,4);
    plot(file.file.t, file.file.uGPSm, "LineWidth",0.5,"color",[0.9290, 0.6940, 0.1250]	); %% inertial!!! Ground speed
    hold on    
    %plot(file.file.t, file.file.XX_k1_k1(4,:)); %% body frame!!!
    plot(file.file.t, file.file.u_inertial_KF+ file.file.XX_k1_k1(16,:),"color","b"); % inertial ground speed
    hold on
    plot(file.file.t,file.file.u_n+ -10, "--","color","r"); %% intertial!!! Ground speed
    title("Velocity u");
    xlabel("Time [s]");
    ylabel("Velocity [m/s]");
    xlim([0 120])
    %legend("Real u","Estimated u", "GPS u");
    
    subplot(3,3,5);
    plot(file.file.t, file.file.vGPSm, "LineWidth",0.5,"color",[0.9290, 0.6940, 0.1250]	); %% inertial!!!   
    hold on 
    %plot(file.file.t, file.file.XX_k1_k1(5,:)); %% body frame!!!
    plot(file.file.t, file.file.v_inertial_KF + file.file.XX_k1_k1(17,:),"color","b");
    hold on
    plot(file.file.t,file.file.v_n+3, "--","color","r"); %% intertial!!!
    title("Velocity v");
    xlabel("Time [s]");
    ylabel("Velocity [m/s]");
    xlim([0 120])
    %legend("Real v","Estimated v", "GPS v");
    
    subplot(3,3,6);
    plot(file.file.t, file.file.wGPSm, "LineWidth",0.5,"color",[0.9290, 0.6940, 0.1250]	); %% inertial!!!
    hold on
    %plot(file.file.t, file.file.XX_k1_k1(6,:)); %% body frame!!!
    plot(file.file.t, file.file.w_inertial_KF + file.file.XX_k1_k1(18,:),"color","b");
    hold on
    plot(file.file.t,file.file.w_n+2, "--","color","r"); %% intertial!!!
    title("Velocity w");
    xlabel("Time [s]");
    ylabel("Velocity [m/s]");
    xlim([0 120])
    %legend("Real w","Estimated w", "GPS w");

    subplot(3,3,7);
    plot(file.file.t, file.file.phiGPSm.*180/pi, "LineWidth",0.5,"color",[0.9290, 0.6940, 0.1250]	);
    hold on 
    plot(file.file.t, file.file.XX_k1_k1(7,:).*180/pi,"color","b");
    hold on
    plot(file.file.t,file.file.phi.*180/pi, "--","color","r");
    title("Angle Phi");
    xlabel("Time [s]");
    ylabel("Angle [deg]");
    xlim([0 120])
    %legend("Real phi","Estimated phi", "GPS phi");

    subplot(3,3,8);
    plot(file.file.t, file.file.thetaGPSm.*180/pi, "LineWidth",0.5,"color",[0.9290, 0.6940, 0.1250]	); 
    hold on 
    plot(file.file.t, file.file.XX_k1_k1(8,:).*180/pi,"color","b");
    hold on
    plot(file.file.t,file.file.theta.*180/pi, "--","color","r");
    title("Angle Theta");
    xlabel("Time [s]");
    ylabel("Angle [deg]");
    xlim([0 120])
    %legend("Real theta","Estimated theta", "GPS theta");

    subplot(3,3,9);
    plot(file.file.t, file.file.psiGPSm.*180/pi, "LineWidth",0.5, "color",[0.9290, 0.6940, 0.1250]	);    
    hold on 
    plot(file.file.t, file.file.XX_k1_k1(9,:).*180/pi,"color","b");   
    hold on
    plot(file.file.t,file.file.psi.*180/pi, "--","color","r");
    title("Angle Psi");
    xlabel("Time [s]");
    ylabel("Angle [deg]");
     xlim([0 120])
    legend("GPS measurement","Estimated state","Real state");
% % 
    %% Estimation error (real)
plotnumber = plotnumber + 1;
figure(plotnumber)
sgtitle(strcat("Real State Estimation Errors of ", file_names{i}));
real_error_deg = file.file.real_error;
real_error_deg(7:9,:) = file.file.real_error(7:9,:).*180/pi;
STD_deg = file.file.STD_x_cor;
STD_deg(7:9,:) = file.file.STD_x_cor(7:9,:)*180/pi;

for j = 1:9
    subplot(3,3,j)
    plot(file.file.t,real_error_deg(j,:));
    hold on 
    plot(file.file.t,STD_deg(j,:));
    hold on 
    plot(file.file.t,-STD_deg(j,:));
    %hold on 
    %plot(file.file.t, file.file.real_error_norm(j,:), "r--");
    title(strcat(state_names{j}," estimation error"));
    xlim([0 120])
    xlabel("Time [s]");
    if (1 <= j) && (j <= 3)
        ylabel("Distance [m]");
        ylim([-1 1]);
    end
    if (4 <= j) && (j <= 6)
        ylabel("Velocity [m/s]");
        ylim([-0.05 0.05])
             
        
    end
    if (7 <= j) && (j <= 9)
        ylabel("Angle [deg]");
        ylim([-0.05 0.05])  
    end
    
    
    %legend('Estimation error', 'Upper error STD', 'Lower error STD');        

end
legend('Estimation error', 'Upper error STD', 'Lower error STD');%, 'Relative error');   
    
    %% data stuff %%
% plotnumber = plotnumber + 1;
% figure(plotnumber)
% plot(file.file.t,file.file.Axm);
% hold on 
% plot(file.file.t, file.file.Aym);
% hold on
% plot(file.file.t, file.file.Azm);
% hold on
% plot(file.file.t, file.file.pm);
% hold on
% plot(file.file.t, file.file.qm);
% hold on
% plot(file.file.t, file.file.rm);
% title(strcat("Input data: ",file_names(i)));
% legend("Axm","Aym","Azm","pm","qm","rm");
% 
% plotnumber = plotnumber + 1;
% figure(plotnumber)
% plot(file.file.t,file.file.Axm);
% hold on 
% plot(file.file.t, file.file.Aym);
% hold on
% plot(file.file.t, file.file.Azm);
% hold on
% plot(file.file.t, file.file.pm);
% hold on
% plot(file.file.t, file.file.qm);
% hold on
% plot(file.file.t, file.file.rm);
% title(strcat("Input data: ",file_names(i)));
% legend("Axm","Aym","Azm","pm","qm","rm");


    %% Innovation %%
% plotnumber = plotnumber + 1;
% figure(plotnumber)
% plot(file.file.t,file.file.innovation);
% title(strcat("Innovation: ",file_names{i}));
% xlabel("Time [s]");
% ylabel("[-]");
% ylim([0 100])
% legend("x_{GPS}","y_{GPS}","z_{GPS}","u_{GPS}","v_{GPS}","w_{GPS}", "phi_{GPS}","theta_{GPS}","psi_{GPS}", "V_{TAS}","alpha","beta");

%% WIND %%
% plotnumber = plotnumber + 1;
% figure(plotnumber)
% plot(file.file.t, file.file.XX_k1_k1(16,:));
% hold on
% plot(file.file.t, file.file.XX_k1_k1(17,:));
% hold on
% plot(file.file.t, file.file.XX_k1_k1(18,:));
% title(strcat("Kalman Filter results of ",file_names{i}));
% xlabel("Time [s]");
% ylabel("Velocity [m/s]");
% legend("Wind X", "Wind Y", "Wind Z");


end

%% Plot Relative Error of all states for all maneuvers in one 3x3 plot %%
% plotnumber = plotnumber + 1;
% figure(plotnumber)
% sgtitle("Relative Estimate Errors for all maneuvers");
% for i = 1:9 %iterate over the states 
%     subplot(3,3,i)
%     for j = 1:length(file_names)
%         dummy = load(strcat("Kalman_data\",file_names{j}));
%         plot(dummy.file.file.t,dummy.file.file.real_error_norm(i,:));
%         hold on
%     end
%     title(strcat(state_names(i), " estimation error"));     
%     
% end
% legend(file_names);

%% BIAS %% 
real_biases = ones(6,length(file.file.t));
real_biases(1,:) = real_biases(1)*0.01;
real_biases(2,:) = real_biases(2)*0.01;
real_biases(3,:) = real_biases(3)*0.01;
real_biases(4,:) = real_biases(4)*0.001;
real_biases(5,:) = real_biases(5)*0.001;
real_biases(6,:) = real_biases(6)*0.001;
bias_names = ["x_{IMU}","y_{IMU}","z_{IMU}","p_{IMU}","q_{IMU}","r_{IMU}"];


plotnumber = plotnumber + 1;
figure(plotnumber)
sgtitle("Bias Estimate Errors of IMU Measurement Data for All Maneuvers");
for i = 1:6 % 6 biases
    average_upper_STD = zeros(1,12001);
    average_lower_STD = zeros(1,12001);
    subplot(2,3,i)
    for j =  1:length(file_names)
        dummy = load(strcat("Kalman_data\",file_names{j}));
        %ytickformat('%.1f')
        if (1 <= i) && (i <= 3)
            plot(dummy.file.file.t,(dummy.file.file.XX_k1_k1(i+9,:)-real_biases(i,:)),"-")%./real_biases(i,:));
            %xlim([0 120])
            hold on
            ylabel("Velocity [m/s]")
            xlabel("Time [s]")
            average_upper_STD = average_upper_STD + dummy.file.file.STD_x_cor(i+9,:);
            average_lower_STD = average_lower_STD - dummy.file.file.STD_x_cor(i+9,:);
            xlim([0 120])
            ylim([-0.008 0.008])
            %legend("upper STD", "lower STD" )
            
        else
            plot(dummy.file.file.t,(dummy.file.file.XX_k1_k1(i+9,:).*180/pi-real_biases(i,:)),"-")%./real_biases(i,:));
            %xlim([0 120])
            hold on
            ylabel("Angular rate [deg/s]")
            xlabel("Time [s]")
            average_upper_STD = average_upper_STD + dummy.file.file.STD_x_cor(i+9,:)*180/pi;
            average_lower_STD = average_lower_STD - dummy.file.file.STD_x_cor(i+9,:)*180/pi;
            xlim([0 120])
            ylim([-0.008 0.008])
            %legend("Upper error STD", "Lower error STD" )
        end
        
    end
    title((bias_names(i)));     
    plot(dummy.file.file.t,average_upper_STD./5, "--", "LineWidth", 1.2)
    hold on
    plot(dummy.file.file.t,average_lower_STD./5,"--","LineWidth", 1.2)
    
end
legend("drdoublet.mat", "dr3211.mat","de3211.mat", "da3211.mat","dadoublet.mat", "Average upper error STD","Average lower error STD");



    