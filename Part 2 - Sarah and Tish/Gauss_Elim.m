%User Input
% x,y,z variables specification
% # of rows
% # of cols
% row input csv values?

clc;
clear all;
close all;
opCount = 0;
row = 3;
col = 3;

n = row

% TODO incorporate user input
A = [1, 2, -1; 2, 1, -2; -3, 1, 1]
ans = [3; 3; -6];
solution = [0; 0; 0];
augA = [A, ans] % tableu form

for j = 1:n-1 % n-1 = num of rows - the 1st row
    % check for division by zero
    if augA(j,j) == 0
        opCount = opCount + 1;
        error('zero pivot encountered');
    end
    % eliminate col j to put zero in each location below diag
    % ex. a(j+1, j), a(j+2, j),...a(n,j)
    for i = j+1:n
        % row multiplier
        multi = augA(i, j) / augA(j, j)
        opCount = opCount + 1;

        % subtract multiplier * the row from 
        for index = 1:col+1
            augA(i, index) = (augA(i, index) - (multi * augA(j, index)))
            opCount = opCount + 1;
            augA(i, index)
        end
    end
end

% Backsolving
for q = n:-1 : 1
    solutions(q) = augA(q, row + 1)
    for u = q+1:n
        opCount = opCount + 1;
        solutions(q) = solutions(q) - augA(q,u)*solutions(u)
    end
    solutions(q) = solutions(q)/augA(q,q)
end

for x = 1:n
    fprintf('x%d = %f\n', x, solutions(x))
end

fprintf('Number of Operations = %d', opCount);

