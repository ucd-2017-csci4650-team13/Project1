% function to calculate forward error, backward error, and error
% magnification of methods used to solve single variable equations
% Pass in original syms function, its derivative, the real root, and the
% approximate root
function calcError(func, deriv, r, xa)
    gPow = feval(symengine, 'degree', func);    % Gets g(x) from highest degree of the equations, func has to be symbolic for degree() to work
    gofr = r^gPow;                              % g(r) will just be r^degree
    magError = abs(gofr/(r*feval(deriv, r)));   % Calculating the magnitude of error from equation on page 49
    forwardErr = abs(r - xa);                   
    backwardErr = abs(subs(func,xa));
    % Remove and return the error magnifcation when integrating code
    fprintf('Real Root = %1.8f \nApproximate Root = %12.8f\nForward Error = %12.8f \nBackward Error = %12.8f\n Error Magnification = %12.8f\n', r, xa, forwardErr, backwardErr, magError);
end


