function [xList, errList, errorFlag, i, opCount, timeStr] = Single_Var_FixedPoint(infxn, x0, r, Tol, maxIterations, handles)
%Tol = 0.00000001;       % Stopping criteria
xList = zeros;          % Have x be a list for creating graphs
errorFlag = false;
errList = zeros;
iCount = zeros;
difference = zeros;
difference(1) = 0;
ticks = 0;
opCount = 0;
f = matlabFunction(infxn);
fprime = matlabFunction(diff(infxn));
% TODO check for convergence
if r ~= Inf
    fprimeofr = fprime(r);
    if abs(fprimeofr) > 1
        divStr = ['Since |g''(r)| = ', num2str(abs(fprimeofr)), ' > 1, FPI may not converge'];
        set(handles.singleVarOutputText, 'string', divStr);
    elseif fprimeofr == 0
        set(handles.singleVarOutputText, 'string', 'Since g''(r) = 0, FPI will be quadratically convergent')
    else
        linConStr = ['Since |g''(r)| < 1, FPI will be linearly convergent with rate ', num2str(abs(fprimeofr))];
        set(handles.singleVarOutputText, 'string', linConStr);
    end
    if r ~= Inf
        errList(1) = abs(xList(1) - r);
    end
end
xList(1) = x0;

% Runs until root approximated or number of iterations reached
tic;
for i = 1:maxIterations
    iCount(i) = i-1;
    xList(i+1) = f(xList(i));           % Get next x from result of current x
    opCount = opCount + 1;
    if r ~= Inf
        errList(i+1) = abs(xList(i) - r);                % Gets forward error of current iteration
    end
    dx = abs(xList(i+1) - xList(i));    % Tracks amount result changed
    opCount = opCount + 1;
    if abs(dx) < Tol
        break;
    end
    difference(i+1) = abs(xList(i+1) - xList(i));
    opCount = opCount + 1;
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
end
timeStr = toc;
