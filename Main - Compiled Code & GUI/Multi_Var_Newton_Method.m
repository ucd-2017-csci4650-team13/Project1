function x = Multi_Var_Newton_Method(vars, x, eqns, number_of_iterations, handles)
%create a Jacobian matrix
DF = jacobian(eqns, vars);
DFStr = char(DF);
set(handles.nonLinearSysOutput, 'string', DFStr);
% disp(DF);
x_values = zeros(1,100);
y_values = zeros(1,100);
t = zeros(1,100);

%begin iteration steps for calculating the solution of the system
for i=1:number_of_iterations
 
    %solve for the solution set, s, to plug into later
    tic;
    a = zeros(length(eqns),1);
    
    for j=1:length(eqns)
        answer = subs(eqns(j), vars, x);
        a(j) = single(answer);
    end %end of solution set loop
    
    
    %solve for the constants in the Jacobian matrix
    %find values of Jacobian with starting point
    sol_matrix = subs(DF, vars, x);
    
    %a = array containing the solution variables [s1, s2, s3...],
    %sol_matrix = solution matrix to solve for the solution set (s1, s2, s3...)
    %solve the linear system of equations this will output the solution set

    
    sol_set = linsolve(sol_matrix, a);
    if(isinf(sol_set))
        warning('The system of equations does not converge')
        return
    end
    sol_set = double(sol_set);
    
    %reshape to a 1D array for easy subtraction
    sol_set = reshape(sol_set, [1, numel(a)]);
    
    
    %solve xk = x(k-1) + s --> xk
    x_values(i) = x(1);
    y_values(i) = x(2);
    prev_x = x;
   
    x = x - sol_set;
    if((round(sum(prev_x - x), 16)) == 0)
        
        figure
        subplot(2,1,1)       % add first plot in 2 x 1 grid
        plot(x_values,y_values)
        title('Convergence or Divergence')
        
        subplot(2,1,2)       % add second plot in 2 x 1 grid
        plot(t)       % plot using + markers
        title('TIme Complexity')
        drawnow()
        return
    end
    diverge = 0;
    if(i <= 10)
        distance(i) = round(sum(prev_x - x));
        if(i == 10)
            if(distance(i) > distance(i-9))
                diverge = 1;
            end
        end
    end
    
    
    %show the solution at the end of each step[[2*u^2 + v^2 + 3*w^2 + 6*w - 4*u + 2],[3*u^2 - 12*u + v^2 + 3*w^2 + 8], [u^2 + v^2 - 2*v + 2*w^2 - 5]]
    t(i) = toc;
end %repeat k times end of iteration loop

if(diverge == 1)
  
    currStr = get(handles.nonLinearSysOutput, 'string');
    newStr = combineString(currStr, 'The solution appears to be diverging');
    set(handles.nonLinearSysOutput, 'string', newStr);
end

figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(x_values,y_values)
title('Convergence or Divergence')

subplot(2,1,2)       % add second plot in 2 x 1 grid
plot(t)       % plot using + markers
title('Time Complexity')