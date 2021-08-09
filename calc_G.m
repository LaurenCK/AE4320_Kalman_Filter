%% Calculates G in x_dot(t) = f(x(t),u(t),t) + G(x(t))*w(t) %%

function [G] = calc_G(x_states)
%x_states = [x, y, z, u, v, w, phi, theta, psi, bias_xr, bias_yr, bias_zr, bias_pr, bias_qr, bias_rr, W_x,W_y,W_z];
u = x_states(4);
v = x_states(5);
w = x_states(6);
phi = x_states(9);
theta = x_states(7);

G = zeros(18,6);
G(4,1) = -1;
G(4,5) = w;
G(4,6) = -v;
G(5,2) = -1;
G(5,4) = -w;
G(5,6) = u;
G(6,3) = -1;
G(6,4) = v;
G(6,5) = -u;
G(7,4) = -1;
G(7,5) = -sin(phi)*tan(theta);
G(7,6) = -cos(phi)*tan(theta);
G(8,5) = -cos(phi);
G(8,6) = sin(phi);
G(9,5) = -sin(phi)/cos(theta);
G(9,6) = -cos(phi)/cos(theta);

end

