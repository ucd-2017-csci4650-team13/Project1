function x1 = Broyden_1_Method(x_i, vars, eqns, A, number_of_iterations)
%evaluate the functions at x
%{(@(u,v)u^2+v^2-1) ,   (@(u,v)(u-1)^2+v^2-1)}

%y1 = cellfun(@(t) t(x{:}), eqns);
y1 = zeros(length(eqns),1);

for j=1:length(eqns)
    answer = subs(eqns(j), vars, x_i);
    y1(j) = single(answer);
end %end of solution set loop
x1 = x_i;
x1 = transpose(x1);

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
 
    y1 = zeros(length(eqns),1);
    for j=1:length(eqns)
        answer = subs(eqns(j), vars, transpose(x1));
        y1(j) = single(answer);
    end %end of solution set loop
    if((round(sum(x - x1), 16)) == 0)
        disp('solution found');
        figure
        plot(t)
        title('Time Complexity Graph')
        xlabel('Iterations')
        ylabel('Time to Compute (s)')
        return
    end
    
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
figure, plot(t)
title('Time vs Iterations')
xlabel('Iterations')
ylabel('Time to Compute (s)')