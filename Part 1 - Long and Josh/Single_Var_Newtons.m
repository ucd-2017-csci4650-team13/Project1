syms x;
Tol = 0.00000001;
%Reading in input as a symbolic variable to use diff
infxn = input('Give an equation in x: ');%,'s');
fxn = matlabFunction(infxn);    % Converting it to a matlab function handle
x = double(input('Enter the initial guess: '));
r = double(input('Enter the real root: '));
infprime = diff(infxn);         % Gets the derivative
fprime = matlabFunction(infprime);  

i = 0;
dx = 1;
pastError = 0;

fprintf('step       x      ei = |xi - r|  ei/(ei-1)^2\n');
fprintf('----   ---------   -----------   ----------\n');

while (dx > Tol || abs(fx) > Tol)   % While difference between xi+1 and xi is above the tolerance
    fx = feval(fxn, x);
    fprimex = feval(fprime, x);
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

backwarderror = abs(fxn(x));
fprintf('Backward error = %12.8f\n', backwarderror);
%TODO Relative Errors and Error Magnification