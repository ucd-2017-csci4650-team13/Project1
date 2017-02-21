%Fixed Point Iteration
% TODO, make into function
syms x;
Tol = 0.00000001;   % Stopping criteria
xList = zeros;          % Have x be a list for creating graphs

infxn = input('Give an equation in x: ');
f = matlabFunction(infxn); % Equation definition

xList(1) = double(input('Enter the initial guess: '));
dx = 1;
i = 1;

while(abs(dx) > Tol)
    xList(i+1) = f(xList(i));           % Get next x from result of current x
    dx = abs(xList(i+1) - xList(i));    % Tracks amount result changed
    i = i + 1;
end

xa = xList(i);
fprintf('Approximate root = %12.8f\n', xa);