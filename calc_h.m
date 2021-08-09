%% Calculates h in z(t) = h(x(t),u(t),t)

function [h_val] = calc_h(x_states)
[x, y, z, u, v, w, phi, theta, psi_, bias_xr, bias_yr, bias_zr, bias_pr, bias_qr, bias_rr, W_x,W_y,W_z] =deal(x_states(1),x_states(2),x_states(3),x_states(4),x_states(5),x_states(6),x_states(7),x_states(8),x_states(9),x_states(10),x_states(11),x_states(12),x_states(13),x_states(14),x_states(15),x_states(16),x_states(17),x_states(18));

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

