
function [xdot_val] = calc_f(t,x_states, u_inputs)
[x, y, z, u, v, w, phi, theta, psi_, bias_xr, bias_yr, bias_zr, bias_pr, bias_qr, bias_rr, W_x,W_y,W_z] =deal(x_states(1),x_states(2),x_states(3),x_states(4),x_states(5),x_states(6),x_states(7),x_states(8),x_states(9),x_states(10),x_states(11),x_states(12),x_states(13),x_states(14),x_states(15),x_states(16),x_states(17),x_states(18));
[Axm, Aym, Azm, pm,qm,rm]= deal(u_inputs(1),u_inputs(2),u_inputs(3),u_inputs(4),u_inputs(5),u_inputs(6));
g = 9.807; % gravitaional constant Earth
xdot_val = zeros(18,1);

xdot_val(1) = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi_)-(v*cos(phi)-w*sin(phi))*sin(psi_)+W_x;
xdot_val(2) = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi_)+(v*cos(phi)-w*sin(phi))*cos(psi_)+W_y;
xdot_val(3) = (-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+W_z);
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



