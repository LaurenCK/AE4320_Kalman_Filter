%% Calculates the Jacobian of the State Transition Equation f(x(t),u(t),t) %%

function [dF] = calc_dF(x_state_values, u_input_values)
persistent dFx;                                 % persistent varibale used to only construct the symbolic Jacobian once
                                                % Then update with new state variable values
                                               
if isempty(dFx)
    [xdot_val, x_states, u_inputs] = def_f();   % get symbolic definitions
    dFx = jacobian(xdot_val, x_states);         % symbolically calculate Jacobian wrt the states
    %display(dFx);
    dFx = matlabFunction(dFx);                  % get numerical function
    %display(dFx);
end
% bias_pr, bias_qr, bias_rr,phi, pm, psi, qm, rm, theta, u, v, w
dF = feval(dFx, x_state_values(13),x_state_values(14),x_state_values(15), x_state_values(7),u_input_values(4),x_state_values(9),u_input_values(5),u_input_values(6),x_state_values(8),x_state_values(4),x_state_values(5),x_state_values(6));
end

