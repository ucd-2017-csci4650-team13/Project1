%User Input
% x,y,z variables specification
% # of rows
% # of cols
% row input csv values?

clc;
clear all;
close all;

row = 3;
col = 3;

n = row

% TODO incorporate user input
A = [1, 2, -1; 2, 1, -2; -3, 1, 1]
ans = [3; 3; 6];

augA = [A, ans]


for j = 1:n-1 % n-1 = num of rows - the 1st row
    % check for division by zero
    if augA(j,j) == 0
        error('zero pivot encountered');
    end
    % eliminate col j to put zero in each location below diag
    % ex. a(j+1, j), a(j+2, j),...a(n,j)
    for i = j+1:n
        % row multiplier
        multi = augA(i, j) / augA(j, j)

        % subtract multiplier * the row from 
        for index = j+1:col+1
            augA(i, index) = augA(i, index) - (multi * augA(j, index))
            augA(i, index)
        end
    end
end


