x_i = input('enter a starting point');
A = input('enter an initial matrix');
eqns = input('enter a system of equations');
number_of_iterations = input('enter the number of iterations');

%evaluate the functions at x
%{(@(u,v)u^2+v^2-1) ,   (@(u,v)(u-1)^2+v^2-1)}
x = num2cell(x_i);
y1 = cellfun(@(t) t(x{:}), eqns);
x1 = x_i;
x1 = transpose(x1);
y1 = transpose(y1);


for i=1:number_of_iterations
    x = x1;
    y = y1;
    x1 = x - A\y;
    x_val = num2cell(x1);
    display(x1);
    y1 = cellfun(@(t) t(x_val{:}), eqns);
    y1 = transpose(y1);
    
    %calculate new matrix A
    deltaY = y1 - y;
    deltaX = x1 - x;
    p1 = A*deltaX;
    p2 = deltaY - p1;
    p3 = p2 * transpose(deltaX);
    p4 = transpose(deltaX)*deltaX;
    A = A + p3/p4;
end




