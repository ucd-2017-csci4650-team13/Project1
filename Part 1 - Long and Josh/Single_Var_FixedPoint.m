%Fixed Point Iteration
% TODO, make into function
syms x;
Tol = 0.00000001;   % Stopping criteria
x = zeros;          % Have x be a list for creating graphs

infxn = input('Give an equation in x: ');
f = matlabFunction(infxn); % Equation definition

x(1) = double(input('Enter the initial guess: '));
dx = 1;
i = 1;

while(abs(dx) > Tol)
    x(i+1) = f(x(i));           % Get next x from result of current x
    dx = abs(x(i+1) - x(i));    % Tracks amount result changed
    i = i + 1;
end

xa = x(i);
fprintf('Approximate root = %12.8f\n', xa);