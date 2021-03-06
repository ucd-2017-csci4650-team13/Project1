%test equations:
%{[@(u,v,w)2*u^2 - 4*u + v^2 + 3*w^2 + 6*w + 2], [@(u,v,w)u^2 +v^2 - 2*v+ 2*w^2 - 5], [@(u,v,w)3*u^2 - 12*u+v^2 + 3*w^2 + 8]}
% initial matrix: [1 0 0; 0 1 0; 0 0 1]
%INPUTS:
% starting point: [1,1]
% initial matrix: [1 0; 0 1]
% system of equations: {(@(u,v)u^2+v^2-1) ,   (@(u,v)(u-1)^2+v^2-1)}
% number of iterations: 2
x_i = input('enter a starting point');
A = input('enter an initial matrix');
eqns = input('enter a system of equations');

%error checking for inputs
if(size(A,1) ~= length(x_i)) 
    disp('the matrix row length must equal the number of variables. Please start again.')
    return
end
if(size(A,2) ~= length(eqns))
    disp('The matrix column length must equal the number of equations. Please start again.')
    return
end
if(length(eqns) ~= length(x_i))
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
x_values = zeros(1,100);
y_values = zeros(1,100);
t = zeros(1,100);
%begin iteration steps
for i=1:number_of_iterations
    x = x1;
    x_values(i) = x(1);
    y_values(i) = x(2);
    y = y1;
    tic;
    x1 = x - A\y;
    disp(x1);
    x_val = num2cell(x1);
    
    y1 = cellfun(@(t) t(x_val{:}), eqns);
    if(round(y1, 10) == 0)
        disp('solution found at x = ')
        disp(vpa(x1,10))
        disp('correct to 10 decimals digits')
        plot(t)
        title('Time Complexity Graph')
        xlabel('Iterations')
        ylabel('Time to Compute (s)')
        return
    end
    y1 = transpose(y1);
    
    %calculate new matrix A
    deltaY = y1 - y;
    deltaX = x1 - x;
    p1 = A*deltaX;
    p2 = deltaY - p1;
    p3 = p2 * transpose(deltaX);
    p4 = transpose(deltaX)*deltaX;
    A = A + p3/p4;
    t(i) = toc;
end
plot(t)



