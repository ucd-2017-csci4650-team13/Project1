%Program 1.1 Bisection Method
%Computes approximate solution of f(x)=0
%Input: function handle f; a,b such that f(a)*f(b)<0,
% and tolerance tol
%Output: Approximate solution xc
function xc=Bisection(infxn,a,b,r,tol)
syms x;
f = matlabFunction(infxn);
if sign(f(a))*sign(f(b)) >= 0
    error('f(a)f(b)<0 not satisfied!') %ceases execution
end
fa=f(a);
fb=f(b);
i = 1;
cList = zeros();
while (b-a)/2>tol
    cList(i)=(a+b)/2;
    fc=f(cList(i));
    if fc == 0 %c is a solution, done
        break
    end
    if sign(fc)*sign(fa)<0 %a and c make the new interval
        b=cList(i);fb=fc;
    else %c and b make the new interval
        a=cList(i);fa=fc;
    end
    i = i + 1;
end
xc=(a+b)/2; %new midpoint is best estimate
calcError(infxn, matlabFunction(diff(infxn)), r, xc)
