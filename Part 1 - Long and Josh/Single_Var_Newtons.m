syms x;
Tol = 0.00000001;
infxn = input('Give an equation in x: ');%,'s');
fxn = matlabFunction(infxn);
x = double(input('Enter the initial guess: '));
infprime = diff(infxn);
fprime = matlabFunction(infprime);
%fx = subs(fxn)
%fprimex = subs(fprime)
%xip1 = x - fx/fprimex
%x = xip1
i = 0;
dx = 1;
fprintf('step       x           dx         f(x)\n');
fprintf('----   ---------   -----------   ----------\n');

while (dx > eps || abs(fx) > eps)
    fx = feval(fxn, x);
    fprimex = feval(fprime, x);
    if(fprimex == 0 || isNaN(fprimex) || abs(fprimex) == Inf)
        error('The derivative of the function at %12.8f is 0, try another initial guess', x);
    end
    xip1 = x - (fx/fprimex);
    fprintf('%3i %12.8f %12.8f %12.8f\n', i, x, dx, xip1);
    dx = abs(x-xip1);
    x = xip1;
    i = i + 1;
end

backwarderror = abs(fxn(x));
fprintf('Backward error = %12.8f\n', backwarderror);