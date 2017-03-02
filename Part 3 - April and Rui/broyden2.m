%INPUTS:
% starting point: [1,1]
% initial matrix: [1 0; 0 1]
% system of equations: {(@(u,v)u^2+v^2-1) ,   (@(u,v)(u-1)^2+v^2-1)}
% number of iterations: 2
x_i = input('enter a starting point');
B = input('enter an initial matrix');
eqns = input('enter a system of equations');

%error checking for inputs
if(size(B,1) < length(x_i)) 
    disp('the matrix row length must equal the number of variables. Please start again.')
    return
end
if(size(B,2) < length(eqns))
    disp('The matrix column length must equal the number of equations. Please start again.')
    return
end
if(length(eqns) < length(x_i))
    disp('the number of equations must equal the number of variables. Please start again.')
    return
end
number_of_iterations = input('enter the number of iterations');

%evaluate the functions at x
%{(@(u,v)u^2+v^2-1) ,   (@(u,v)(u-1)^2+v^2-1)}
x = num2cell(x_i);
y1 = cellfun(@(t) t(x{:}), eqns);
x1 = x_i;
x1 = transpose(x1);
y1 = transpose(y1);

%begin iteration steps
for i=1:number_of_iterations
    x = x1;
    y = y1;
    x1 = x - B*y
    x_val = num2cell(x1);
    display(x1);
    y1 = cellfun(@(t) t(x_val{:}), eqns);
    if(round(y1, 10) == 0)
        disp('solution found at x = ')
        disp(vpa(x1,10))
        disp('correct to 10 decimals digits')
        return
    end
    y1 = transpose(y1);
    
    %calculate new matrix B
    deltaY = y1 - y;
    deltaX = x1 - x;
    p1 = B*deltaY;
    p2 = deltaX - p1;
    p3 = p2 * (transpose(deltaX)*B);
    p4 = transpose(deltaX)*B*deltaY;
    B = B + p3/p4;
end




