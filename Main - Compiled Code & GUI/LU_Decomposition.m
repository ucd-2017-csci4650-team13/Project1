% This function performs the LU factorization of a symmetric matrix.
% The function checks for singularity by checking for a zero pivot.
% The function also does a partial pivoting protocol to prevent swamping.
% This function also creates a permutation matrix to help with backsolving.
%   Input: Symmetric matrix and solution vector
%   Output: Lower Triangular factor
%           Upper Triangular factor
%           Permutation matrix
%           Solution Vector x

% input from user for now. may need to update.
function [solution, errFlag] = LU_Decomposition(augA, handles)
%augA = [1 0 1; 0 1 1,
rows = size(augA,1);
columns = rows + 1;
A = augA;
A(:,columns) = [];
%A = input('enter an initial matrix: \n');
b = augA(:,columns);
opiCount = 0;
errFlag = false;

% size of n x n
[~, n] = size(A);

% init U to A
U = A;

% init L to identity
L = eye(n);

% identity matrix for pivoting
P = eye(n);

tic;
for k = 1:n
    % create matrix pivot
    %  maximum absolute value from the k thru n rows, k column
    [pivot i] = max(abs(U(k:n,k)));
    i = i + k-1;
    opiCount = opiCount + 1;
    % if i != k
    if i ~= k
        % interchange rows i and k in U
        % U(k, :) means (row k, all columns in row k)
        opiCount = opiCount + 3;
        temp = U(k, :);
        U(k, :) = U(i, :);
        U(i, :) = temp;
        % interchange rows i and k in P
        opiCount = opiCount + 3;
        temp = P(k, :);
        P(k, :) = P(i, :);
        P(i, :) = temp;
        
        % Lower triangular array
        if k >= 2
            % (k , 1:k-1) means entire row k and columns 1 thru k-1
            temp = L(k, 1:k-1);
            L(k, 1:k-1) = L(i, 1:k-1);
            L(i, 1:k-1) = temp;
            opiCount = opiCount + 3;
        end %endif
    end %end outer if ~=
    % get factors L and U
    for j = k+1:n
        % set multiplier in L
        L(j, k) = U(j, k) / U(k, k);
        opiCount = opiCount + 1;
        % for all in row j = for all in row j - multiplier*for all in row k
        U(j, :) = U(j, :) - L(j, k) * U(k, :);
        opiCount = opiCount + 1;
    end % end gauss elim
end

%solve for c
c = [];
inv(L)
Pb = P*b;
c = L\Pb;
opiCount = opiCount + n^2;

%solve for x
solution = [];
solution = U\c;
opiCount = opiCount + n^2;
solution = solution';
time = toc;

opsString = ['Number of Operations = ', num2str(opiCount)];
timeString = ['Elapsed Time = ', num2str(time), ' seconds'];
newString = combineString(opsString, timeString);
set(handles.lSysOutputText, 'string', newString);