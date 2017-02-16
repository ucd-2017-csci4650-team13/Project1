%TODO, make into function
syms x;
Tol = 0.00000001;
%Reading in input as a symbolic variable to use diff
infxn = input('Give an equation in x: ');
x = double(input('Enter the initial guess: '));
r = double(input('Enter the real root: '));
f = matlabFunction(infxn)  ;            % Equation definition
fprime = matlabFunction(diff(infxn));   % First-order derivative of f
%fprime = diff(fxn);                    % Gets the derivative 

i = 0;                                  % Tracks the iteration
dx = 1;
fofx = 1;
pastError = 0;

fprintf('step       x      ei = |xi - r|  ei/(ei-1)^2\n');
fprintf('----   ---------   -----------   ----------\n');

while (dx > Tol || abs(fofx) > Tol)   % While difference between xi+1 and xi is above the tolerance
    fofx = f(x);
    fpofx = fprime(x);
    if(fpofx == 0 || abs(fpofx) == Inf)
        error('The derivative of the function at %12.8f is 0, try another initial guess', x);
    end
    
    xip1 = x - fofx/fpofx;              % Gets the next value of x
    ei = abs(x - r);                    % Gets forward error of current iteration
    iterativeError = ei/(pastError)^2;  % Error in relation to previous iteration
    
    if iterativeError == Inf || iterativeError > 1  % Handles case where result is invalid
        iterativeError = NaN;
    end
    
    fprintf('%3i %12.8f %12.8f   %12.8f %12.8f\n', i, x, ei, iterativeError, xip1);
    
    % Prep for next iteration of loop
    pastError = ei;
    x = xip1;
    i = i + 1;
end

fprintf('root = %12.8f\n', x);
calcError(infxn, fprime, r, x);
fprintf('root %12.8f has m = %12.8f\n', r, getRootMultiplicity(infxn, r));