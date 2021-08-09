%% Calculate and plot the change in KF results due to 10x higher variance V, alpha, beta %%
clear("all")
close("all")
file_names = ["drdoublet.mat", "dr3211.mat","de3211.mat", "da3211.mat","drdoublet.mat"];
state_names = {'x','y','z','u','v','w','Phi','Theta','Psi'};
bias_names = ["x_{IMU}","y_{IMU}","z_{IMU}","p_{IMU}","q_{IMU}","r_{IMU}"];
real_biases = ones(6,12001);
real_biases(1,:) = real_biases(1)*0.01;
real_biases(2,:) = real_biases(2)*0.01;
real_biases(3,:) = real_biases(3)*0.01;
real_biases(4,:) = real_biases(4)*0.001*pi/180;
real_biases(5,:) = real_biases(5)*0.001*pi/180;
real_biases(6,:) = real_biases(6)*0.001*pi/180;


%data_generation(0);
%Kalman(0.1,0.1*pi/180,0.1*pi/180,0);
%data_generation(1);
%Kalman(0.1,0.1*pi/180,0.1*pi/180,1); %save in Kalman_data2


plotnumber = 0;
%% Make plot of state estimation error %%
for i = 1:1%length(file_names)
    file1 = load(strcat("Kalman_data\",file_names{i}));
    file2 = load(strcat("Kalman_data2\",file_names{i}));
    % states
    plotnumber = plotnumber + 1;
    figure(plotnumber);
    sgtitle(strcat("State Estimation Errors for ", file_names(i)));

    for j = 1:9 % 9 states
        subplot(3,3,j);
        plot(file1.file.file.t,file1.file.file.real_error(j,:));
        hold on 
        plot(file2.file.file.t,file2.file.file.real_error(j,:));
        hold on
        plot(file1.file.file.t, file1.file.file.STD_x_cor(j,:)), "--";
        hold on
        plot(file1.file.file.t, -1*file1.file.file.STD_x_cor(j,:), "--");
        title(strcat(state_names{j}));
        xlabel("Time [s]")
        grid on
        if (1 <= j) && (j<=3)
            ylabel("Distance [m]")
            ylim([-0.5,0.5])
        end
        
        if (4 <= j) && (j<=6)
            ylabel("Velocity [m/s]")
            ylim([-0.1,0.1])
            
        end
        if (7 <= j) && (j<=9)
            ylabel("Angle [rad]")
            ylim([-0.001,0.001])
        end
        
    end
    legend("Run 1", "Run 2", "Average upper error STD","Average lower error STD");

    % bias
    plotnumber = plotnumber + 1;
    figure(plotnumber);
    sgtitle(strcat("IMU Bias Estimation Errors for ", file_names(i)));

    for j = 1:6 % 6 biases
        subplot(2,3,j);
        plot(file1.file.file.t,(file1.file.file.XX_k1_k1(j+9,:)-real_biases(j,:)));
        hold on 
        plot(file1.file.file.t,(file2.file.file.XX_k1_k1(j+9,:)-real_biases(j,:)));
        hold on
        plot(file1.file.file.t,(file1.file.file.STD_x_cor(j+9,:)),"--");
        hold on
        plot(file1.file.file.t,(-1*file1.file.file.STD_x_cor(j+9,:)),"--");
        title(strcat(bias_names{j}));
        xlabel("Time [s]")
        if (1<=j)&&(j<=3)
            ylim([-0.01 0.01])
            ylabel("Velocity [m/s]")
        else
            ylim([-10e-5 10e-5]) 
            ylabel("Angular rate [rad/s]")
        end
        grid on
    end
    legend("Run 1", "Run 2","Average upper error STD","Average lower error STD");
    


end


