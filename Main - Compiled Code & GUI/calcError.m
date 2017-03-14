% function to calculate forward error, backward error, and error
% magnification of methods used to solve single variable equations
% Pass in original syms function, its derivative, the real root, and the
% approximate root
function calcError(func, deriv, r, xa, handles, i)
derivErrorFlag = false;
% might have to remove error magnification or find alternative for gpow

backwardErr = double(abs(subs(func,xa)));
approxStr = 'Approximate Root = ';
if r ~= Inf
    forwardErr = abs(r - xa);
    rootM = getRootMultiplicity(func, r);
    realString = 'Real Root = ';
    mString = ' has multiplicity = ';
    fString = 'Forward Error = ';
    forwardStr = sprintf('%s%s%s%s\n%s%s', realString, num2str(r, '%20.10f'), mString, num2str(rootM, '%20.10f'), fString, num2str(forwardErr, '%20.10f'));
    try
        gPow = feval(symengine, 'degree', func);    % Gets g(x) from highest degree of the equations, func has to be symbolic for degree() to work
        gofr = r^gPow;                              % g(r) will just be r^degree
        magError = double(abs(gofr/(r*feval(deriv, r))));   % Calculating the magnitude of error from equation on page 49
    catch
        derivErrorFlag = true;
    end
end

iterationsString = ['Number of Iterations = ', num2str(i)];
bString = 'Backward Error = ';
currString = get(handles.singleVarOutputText, 'string');
newString = sprintf('%s%s\n%s%s%s%s\n%s%s\n%s%s\n', approxStr, num2str(xa, '%20.10f'), bString, num2str(backwardErr, '%20.10f'));
statusString = combineString(currString, iterationsString);
statusString = combineString(statusString, newString);
if r ~= Inf
    statusString = combineString(statusString, forwardStr);
end
if derivErrorFlag == false && r ~= Inf
    digits = floor(log10(magError*10));
    errorMagStr = ['The function has an error magnification of: ', num2str(magError), ', ', num2str(digits), ' digits of accuracy are lost.'];
    statusString = combineString(statusString, errorMagStr);
end
set(handles.singleVarOutputText, 'string', statusString);
% Remove and return the error magnifcation when integrating code