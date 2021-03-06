function multiplicity = getRootMultiplicity(fnx, x)
multiplicity = 0;
while(subs(fnx, x) == 0)        % Runs when fn(0) = 0
    multiplicity = multiplicity + 1;
    fnx = diff(fnx);            % Get the next derivative
end