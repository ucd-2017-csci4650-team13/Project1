function [cList, errList, errorFlag, i, opCount, timeStr] = Bisection(infxn,a,b,r,tol, maxIterations, handles)
f = matlabFunction(infxn);
errorFlag = false;
iCount = zeros;
cList = zeros;
opCount = 0;

errList = zeros;
i = 0;
tic;
if sign(f(a))*sign(f(b)) >= 0
    errorString = 'f(a)f(b)<0, Intermediate Value Theorem not satisfied!'; %ceases exe  cution
    set(handles.singleVarOutputText, 'string', errorString);
    errorFlag = true;
else
    cList = zeros;
    fa=f(a);
    
    %fb=f(b);
    for i = 1:maxIterations
        iCount(i) = i-1;
        cList(i)=(a+b)/2;
        opCount = opCount + 1;
        if r ~= Inf
            errList(i) = abs(cList(i) - r);                % Gets forward error of current iteration
        end
        
        fc=f(cList(i));
        if fc == 0 || (b-a)/2<tol               %c is a solution or below the tolerance, done
            break
        end
        if sign(fc)*sign(fa)<0 %a and c make the new interval
            b=cList(i);%fb=fc;
        else %c and b make the new interval
            a=cList(i);fa=fc;
        end
        opCount = opCount + 1;
    end
    cList(i)=(a+b)/2; %new midpoint is best estimate
    opCount = opCount + 1;
end
i = i - 1; % To account for extra iteration.
timeStr = toc;