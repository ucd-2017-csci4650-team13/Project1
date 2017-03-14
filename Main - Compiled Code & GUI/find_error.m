function [forward, backward]=find_error(xa, x, eqns,vars)
%forward error

forward = norm((x - transpose(xa)),inf);

%backward error
y1 = zeros(1,length(eqns));
for j=1:length(eqns)
    answer = subs(eqns(j), vars, x);
    y1(j) = single(answer);
end %end of solution set loop
backward = norm(y1, inf);