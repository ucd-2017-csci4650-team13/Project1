%TODO, make into function
syms x;
Tol = 0.00000001;
%Reading in input as a symbolic variable to use diff
fxn = input('Give an equation in x: ');
x = double(input('Enter the initial guess: '));
r = double(input('Enter the real root: '));
fprime = diff(fxn);                 % Gets the derivative 

i = 0;                              % Tracks the iteration
dx = 1;
pastError = 0;

fprintf('step       x      ei = |xi - r|  ei/(ei-1)^2\n');
fprintf('----   ---------   -----------   ----------\n');

while (dx > Tol || abs(fx) > Tol)   % While difference between xi+1 and xi is above the tolerance
    fx = subs(fxn, x);
    fprimex = subs(fprime, x);
    if(fprimex == 0 || abs(fprimex) == Inf)
        error('The derivative of the function at %12.8f is 0, try another initial guess', x);
    end
    xip1 = x - (fx/fprimex);
    ei = abs(x - r);
    iterativeError = ei/(pastError)^2;
    if iterativeError == Inf || iterativeError > 1
        iterativeError = NaN;
    end
    fprintf('%3i %12.8f %12.8f   %12.8f %12.8f\n', i, x, ei, iterativeError, xip1);
    % Prep for next iteration of loop
    pastError = ei;
    dx = abs(x-xip1);
    x = xip1;
    i = i + 1;
end

calcError(fxn, fprime, r, x);
%TODO Relative Errors and Error Magnification