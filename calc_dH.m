%% Calculates the Jacobian of the observation equation Hx = d/dx( h(x(t),u(t),t) ) %%

function [dH] = calc_dH(x_state_values)
persistent dHx
if isempty(dHx) % Use of persistent variable in order to only calculate the Jacobian once 
                % Then only update state input values
 
    %% Calculate Jacobian %%
    [h_val, x_states] = def_h();    % Get symbolic definitions
    dHx = jacobian(h_val,x_states); % Caculate Jacobian symbollically
    %display(dHx);
    dHx = matlabFunction(dHx);      % Make numerical function
    %display(dHx);
    
end

%input = num2cell(x_state_values);   % make cell to be able to select variables individually using {:}
%dH = feval(dHc, input{:}); % gives error that input has too may arguments
% because dHc only has phi, theta, psi, u,v,w as variables...

% manually select states phi, psi,theta, u, v, w as states:
dH = feval(dHx, x_state_values(7), x_state_values(9),x_state_values(8),x_state_values(4),x_state_values(5),x_state_values(6));          % calculates the dHx function on the input state values         
                                    % feval: evaluates a function using its name or its handle, and using the input arguments
                                    % input{:} takes all variable seperately as inputs

end

