%to test function I used the following inputs (this from the text book page 132):
% variables -> {'u', 'v'}
% intial starting point -> [1,2]
% array of equations -> [(v + (-u^3)), (u^2+v^2+(-1))]
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
v = input('enter the variables used for the systems in {''v1'', ''v2''} format');

%convert the cell array to a sym array and each element to sym class
vars = cell2sym(v);
for i=1:numel(v)
    syms(v(i))
end

%Ask user for:
% 1. x = initial starting point,
% 2. r = the array of equations and the
% 3. number of iterations to run through
x = input('enter the initial starting point in [x0, x1] format');

r = input('enter in the array of equations in [(eq1), (eq2), (eq3)] format');

number_of_iterations = input('enter the number of iterations');

%create a Jacobian matrix
DF = jacobian(r, vars);

%begin iteration steps for calculating the solution of the system
for i=1:number_of_iterations
    
    %solve for the solution set, s, to plug into later
    a = zeros(numel(r),1);
    for j=1:numel(r)
        answer = subs(r(j), vars, x);
        a(j) = answer;
    end %end of solution set loop


    %solve for the constants in the Jacobian matrix
    %find values of Jacobian with starting point
    sol_matrix = subs(DF, vars, x);

    %a = array containing the solution variables [s1, s2, s3...],
    %sol_matrix = solution matrix to solve for the solution set (s1, s2, s3...)
    %solve the linear system of equations this will output the solution set
    sol_set = linsolve(sol_matrix, a);
    
    %reshape to a 1D array for easy subtraction
    sol_set = reshape(sol_set, [1, numel(a)]);


    %solve xk = x(k-1) + s --> xk
    x = x - sol_set;
    
    %show the solution at the end of each step
    display(vpa(double(x)))
end %repeat k times end of iteration loop


