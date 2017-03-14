function x1 = Broyden_2_Method(x_i, vars, eqns, B, number_of_iterations)
%evaluate the functions at x

disp(x_i)
y1 = zeros(length(eqns),1);
for j=1:length(eqns)
    answer = subs(eqns(j), vars, x_i);
    y1(j) = single(answer);
end %end of solution set loop

x1 = x_i;
x1 = transpose(x1);

t = zeros(1,100);
%begin iteration steps
for i=1:number_of_iterations
    x = x1;
    y = y1;
    tic;
    x1 = x - B*y;
    
    %display(x1);
    y1 = zeros(length(eqns),1);
    for j=1:length(eqns)
        answer = subs(eqns(j), vars, transpose(x1));
        y1(j) = single(answer);
    end %end of solution set loop
    if((round(sum(x - x1), 16)) == 0)
        %disp('solution found')
        figure
        plot(t)
        title('Time Complexity Graph')
        xlabel('Iterations')
        ylabel('Time to Compute (s)')
        drawnow()
        return
    end
  
    %calculate new matrix B
    deltaY = y1 - y;
    deltaX = x1 - x;
    p1 = B*deltaY;
    p2 = deltaX - p1;
    p3 = p2 * (transpose(deltaX)*B);
    p4 = transpose(deltaX)*B*deltaY;
    B = B + p3/p4;
    t(i) = toc;
end
figure, plot(t)
title('Time vs Iterations')
xlabel('Iterations')
ylabel('Time to Compute (s)')