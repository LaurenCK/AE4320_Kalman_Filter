%% Calculates f in x_dot(t) = f(x(t),u(t),t) + G(x(t))*w(t) %%

function [xdot_val, x_states, u_inputs] = def_f()
syms x y z u v w phi theta psi_ bias_xr bias_yr bias_zr bias_pr bias_qr bias_rr W_x W_y W_z Axm Aym Azm pm qm rm
x_states = [x, y, z, u, v, w, phi, theta, psi_, bias_xr, bias_yr, bias_zr, bias_pr, bias_qr, bias_rr, W_x,W_y,W_z];
u_inputs = [Axm, Aym, Azm, pm,qm,rm];
g = 9.807; % gravitaional constant Earth
xdot_val = sym(zeros(18,1));

xdot_val(1) = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi_)-(v*cos(phi)-w*sin(phi))*sin(psi_)+W_x;
xdot_val(2) = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi_)+(v*cos(phi)-w*sin(phi))*cos(psi_)+W_y;
xdot_val(3) = -u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+W_z;
xdot_val(4) = Axm - bias_xr-g*sin(theta)+(rm-bias_rr)*v-(qm - bias_qr)*w;
xdot_val(5) = Aym - bias_yr +g*cos(theta)*sin(phi) + (pm - bias_pr)*w - (rm - bias_rr)*u;
xdot_val(6) = Azm - bias_zr + g*cos(theta)*cos(phi)+(qm - bias_qr)*u-(pm-bias_pr)*v;
xdot_val(7) = pm - bias_pr +(qm-bias_qr)*sin(phi)*tan(theta)+(rm-bias_rr)*cos(phi)*tan(theta);
xdot_val(8) = (qm - bias_qr)*cos(phi)-(rm-bias_rr)*sin(phi);
xdot_val(9) = (qm-bias_qr)*sin(phi)/cos(theta)+(rm-bias_rr)*cos(phi)/cos(theta);
xdot_val(10) = 0;
xdot_val(11) = 0;
xdot_val(12) = 0;
xdot_val(13) = 0;
xdot_val(14) = 0;
xdot_val(15) = 0;
xdot_val(16) = 0;
xdot_val(17) = 0;
xdot_val(18) = 0;
end

