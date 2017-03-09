%to test function I used the following inputs (this from the text book page 132):
% variables -> {'u', 'v'} or {'u', 'v', 'w'}
% intial starting point -> [1,2] or [1,1,1]
% array of equations -> [[v + (-u^3)], [u^2+v^2+(-1)]] or 
% [[2*u^2 + v^2 + 3*w^2 + 6*w - 4*u + 2],[3*u^2 - 12*u + v^2 + 3*w^2 + 8], [u^2 + v^2 - 2*v + 2*w^2 - 5]]
% number of iterations -> 7
% OUTPUT:
%  
% [ 1.0, 1.0]
%  
% [ 0.875, 0.625]
%  
% [ 0.82903634826711749788672865595943, 0.5643491124260355029585798816568]
%  
% [ 0.82604010817065232075151470780838, 0.56361977350284431231841608678224]
%  
% [ 0.82603135773240976558184911482385, 0.56362416213162991329710393983987]
%  
% [ 0.82603135765418700398043938548653, 0.56362416216125854617757795494981]
%  
% [ 0.82603135765418700398043938548653, 0.56362416216125854617757795494981]


%define variables based on number put in
variables = input('enter the variables used for the systems in {''v1'', ''v2''} format');
disp(numel(variables));
%convert the cell array to a sym array and each element to sym class
vars = cell2sym(variables);
for i=1:length(variables)
     syms(variables(i))
end

%Ask user for:
% 1. x = initial starting point
% 2. eqns = the array of equations and the
% 3. number of iterations to run through
x = input('enter the initial starting point in [x0, x1] format');



eqns = input('enter in the array of equations in [(eq1), (eq2), (eq3)] format');

% if(length(v) ~= length(eqns))
%     warning('There must be the same amount of variables as equations')
%     return
% end

number_of_iterations = input('enter the number of iterations');

%create a Jacobian matrix
DF = jacobian(eqns, vars);
x_values = zeros(1,100);
y_values = zeros(1,100);
t = zeros(1,100);
%begin iteration steps for calculating the solution of the system
for i=1:number_of_iterations
    disp('iteration: ');
    disp(i);
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
    %{'a', 'b'}
    %[1,1]

    %[(a+b), (a+(-b)), (a+b^2)]
    %
    
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
    disp(vpa(x,10))
    x = x - sol_set;
    if((round(sum(prev_x - x), 16)) == 0)
        disp('solution found at x = ');
        disp(vpa(x))
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
    disp('the solution appears to be diverging')
    
    
end
% plot(x_values, y_values)
% drawnow()
% plot(t)
% drawnow()
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(x_values,y_values)
title('Convergence or Divergence')

subplot(2,1,2)       % add second plot in 2 x 1 grid
plot(t)       % plot using + markers
title('Time Complexity')