% Naive Gaussian Elimination
% x,y,z variables specification
% # of rows
% # of cols
% row input csv values?
function [solutions, errFlag] = Gauss_Elim(augA, handles)
errFlag = false;
solutions = zeros;
n = size(augA, 1);
row = size(augA, 1);
col = row + 1;
opCount = 0;
A = augA;
A(:,col) = [];

% CONDITION NUMBER
% Conditioning is a property of the matrix
% Error Magnification factors of the magnitude cond(matrix) are possible
% Matlab default precision is double
norminf = norm(A, inf);
norminf_inv = norm(inv(A), inf);

cond_num = norminf * norminf_inv;

iso_exp = floor(log10(cond_num*10));

condStr1 = ['Error Magnification factors of the magnitude ', num2str(iso_exp), ' are possible'];
condStr2 = 'Since Matlab defaults to double precision this means that';
condStr3 = ['16 - ', num2str(iso_exp),' = ', num2str(16-iso_exp), ' correct digits in the solution.'];

condStr = sprintf('%s\n%s\n%s\n', condStr1, condStr2, condStr3);

set(handles.lSysOutputText, 'string', condStr);

tic;

for j = 1:n-1 % n-1 = num of rows - the 1st row
    % check for division by zero
    if augA(j,j) == 0
        errorStr = 'Zero pivot encountered, stopping.';
        set(handles.lSysOutputText, 'string', errorStr);
        errFlag = true;
        toc;
        break;
    end
    % eliminate col j to put zero in each location below diag
    % ex. a(j+1, j), a(j+2, j),...a(n,j)
    for i = j+1:n
        % row multiplier
        multi = augA(i, j) / augA(j, j);
        opCount = opCount + 1;
        
        % subtract multiplier * the row from
        for index = 1:col
            augA(i, index) = (augA(i, index) - (multi * augA(j, index)));
            opCount = opCount + 2;
        end
    end
end

if errFlag == false
    % Backsolving
    for q = n:-1 : 1
        solutions(q) = augA(q, row + 1);
        for u = q+1:n
            opCount = opCount + 2;
            solutions(q) = solutions(q) - augA(q,u)*solutions(u);
        end
        solutions(q) = solutions(q)/augA(q,q);
        opCount = opCount + 1;
    end
    
    time = toc;
    opsString = ['Number of Operations = ', num2str(opCount)];
    timeString = ['Elapsed Time = ', num2str(time), ' seconds'];
    newString = combineString(condStr, opsString);
    newString = combineString(newString, timeString);
    set(handles.lSysOutputText, 'string', newString);
end