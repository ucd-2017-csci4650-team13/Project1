% Function to implement Single Variable Newton's Method
% Receives:
% infxn = Function in x
% x0 = initial guess
% r = real root
% iterations = number of iterations
% Calculates error and root multiplicity at the end
% Returns a list of x values

function xList = Single_Var_Newtons(infxn, x0, r, iterations)
Tol = 0.00000001;
xList = zeros();                        % List of x values calculated for graphing and tables
xList(1) = x0;                          % Set the first element of the list to the initial guess
iterativeErrorList = zeros();           % List of ei

f = matlabFunction(infxn);              % Convert the symbolic function to a function handle
fprime = matlabFunction(diff(infxn));   % First derivative of f

% For analysis
pastError = 0;
dx = 1;

if(fprime(r) == 0)
    fprintf('f''(r) = 0 so Newton''s Method will be linearly convergent\n');
else
    fprintf('The method will converge quadratically.\n');
end
for i = 1:iterations
    fofx = f(xList(i));
    fpofx = fprime(xList(i));
    
    if(fpofx == 0 || abs(fpofx) == Inf)
        error('The derivative of the function at %12.8f is 0, try another initial guess', xList);
    end
    
    xList(i+1) = xList(i) - fofx/fpofx ;             % Gets the next value of x
    ei = abs(xList(i) - r);                    % Gets forward error of current iteration
    iterativeErrorList(i) = ei/(pastError)^2;  % Error in relation to previous iteration
    
    % Checks if the difference in x values has converged
    if (dx <= Tol || abs(fofx) <= Tol)
        break;
    end
    % Checks if error still exists
    if iterativeErrorList(i) == Inf || iterativeErrorList(i) > 1  % Handles case where result is invalid
        iterativeErrorList(i) = NaN;
    end
        
    % Prep for next iteration of loop
    pastError = ei;
end

fprintf('root = %12.8f\n', xList(i));
calcError(infxn, fprime, r, xList(i));
fprintf('root %12.8f has m = %8f\n', r, getRootMultiplicity(infxn, r));
end