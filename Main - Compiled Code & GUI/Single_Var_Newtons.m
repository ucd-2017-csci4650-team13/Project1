function [xList, errList, errorFlag, i] = Single_Var_Newtons(infxn, x0, r, Tol, maxIterations, handles)
xList = zeros();                        % List of x values calculated for graphing and tables
errList = zeros();
xList(1) = x0;                          % Set the first element of the list to the initial guess
errList(1) = abs(xList(1) - r);         % List of ei
errorFlag = false;
f = matlabFunction(infxn);              % Convert the symbolic function to a function handle
fprime = matlabFunction(diff(infxn));   % First derivative of f
iCount = zeros();
difference = zeros();
ticks = 0;
dx = 1;

if r ~= Inf
    if(fprime(r) == 0)
        fprintf('f''(r) = 0 so Newton''s Method will be linearly convergent\n');
        linearString = '\n f''(r) = 0 so Newton''s Method will be linearly convergent.';
        set(handles.singleVarOutputText, 'string', linearString);
    else
        fprintf('The method will converge quadratically.');
        quadConString = 'The method will converge quadratically.';
        set(handles.singleVarOutputText, 'string', quadConString);
    end
end

for i = 1:maxIterations
    iCount(i) = i-1;
    fofx = f(xList(i));
    fpofx = fprime(xList(i));
    
    if(fpofx == 0 || abs(fpofx) == Inf)
        errorString = ['The derivative of the function at ', num2str(xList(i)), ' is 0, try another initial guess'];
        set(handles.singleVarOutputText, 'string', errorString);
        errorFlag = true;
        break;
    end
    
    xList(i+1) = xList(i) - fofx/fpofx ;             % Gets the next value of x
    difference(i+1) = abs(xList(i+1) - xList(i));
    if (difference(i) > 2 && difference(i) < difference(i + 1))
        ticks = ticks + 1;
    end
    
    if (ticks > 3)
        currStr = get(handles.singleVarOutputText, 'string');
        divStr = 'The answer is diverging.';
        newStr = combineString(currStr, divStr);
        set(handles.singleVarOutputText, 'string', newStr);
        errorFlag = true;
        break;
    end
    
    if r ~= Inf
        errList(i+1) = abs(xList(i+1) - r);                % Gets forward error of current iteration
    end
    % Checks if the difference in x values has converged
    if (dx <= Tol || abs(fofx) <= Tol)
        break;
    end
end