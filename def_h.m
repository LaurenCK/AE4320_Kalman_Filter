%% Calculates h in z(t) = h(x(t),u(t),t)

function [h_val, x_states] = def_h()
syms x y z u v w phi theta psi_ bias_xr bias_yr bias_zr bias_pr bias_qr bias_rr W_x W_y W_z Axm Aym Azm pm qm rm
x_states = [x, y, z, u, v, w, phi, theta, psi_, bias_xr, bias_yr, bias_zr, bias_pr, bias_qr, bias_rr, W_x,W_y,W_z];

x_GPS = x;
y_GPS = y;
z_GPS = z;
u_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi_)-(v*cos(phi)-w*sin(phi))*sin(psi_)+W_x;
v_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi_)+(v*cos(phi)-w*sin(phi))*cos(psi_)+W_y;
w_GPS = -u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+W_z;
phi_GPS = phi;
theta_GPS = theta;
psi_GPS = psi_;
V_TAS = sqrt(u*u +v*v+w*w);
alpha = atan(w/u);
beta = atan(v/sqrt(u*u+w*w));

h_val = [x_GPS,y_GPS,z_GPS,u_GPS,v_GPS,w_GPS,phi_GPS,theta_GPS,psi_GPS,V_TAS,alpha,beta];
end

