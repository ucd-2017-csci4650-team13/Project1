function e=find_error(r, xa, eqns)
%forward error
forward = abs(r - xa);
formatSpec = 'Forward error is %4.2f';
fprintf(formatSpec,forward);

%backward error
x = num2cell(xa);
backward = cellfun(@(t) t(x{:}), eqns);
backward = 
formatSpec = 'Backward error is %4.2f';
fprintf(formatSpec,backward);

%magnification error

end