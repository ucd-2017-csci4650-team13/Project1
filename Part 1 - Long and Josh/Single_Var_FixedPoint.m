% Fixed Point Iteration
% Receives:
% infxn = function in x
% x0 = initial guess
% iterations = number of iterations
% Returns list of calculated x values

function xList = Single_Var_FixedPoint(infxn, x0, iterations) 
Tol = 0.00000001;       % Stopping criteria
xList = zeros;          % Have x be a list for creating graphs

f = matlabFunction(infxn); 

xList(1) = x0;

% Runs until root approximated or number of iterations reached
for i = 1:iterations
    xList(i+1) = f(xList(i));           % Get next x from result of current x
    dx = abs(xList(i+1) - xList(i));    % Tracks amount result changed
    if abs(dx) < Tol
        break;
    end
end

xa = xList(i);
fprintf('Approximate root = %8f\n', xa);
end