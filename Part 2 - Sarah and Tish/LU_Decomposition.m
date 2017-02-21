clc;
clear all;
close all;

% input from user for now. may need to update.
A = input('enter an initial matrix');
b = input('enter an initial vector');

% % % matrix A
% A=[-1, 0, 1;
%    2, 1, 1;
%   -1, 2, 0]
 %vector b
% b = [-2; 17; 3]

% size of n x n
[n, n] = size(A);

% init U to A
U = A;

% init L to identity
L = eye(n);

% identity matrix for pivoting
P = eye(n);

for k = 1:n
    % create matrix pivot
    %  maximum absolute value from the k thru n rows, k column
    [pivot i] = max(abs(U(k:n,k)));
    i = i + k-1;
    
    % if i != k
    if i ~= k 
        % interchange rows i and k in U
        % U(k, :) means (row k, all columns in row k)
        temp = U(k, :);
        U(k, :) = U(i, :);
        U(i, :) = temp;
        % interchange rows i and k in P
        temp = P(k, :);
        P(k, :) = P(i, :);
        P(i, :) = temp;
        
        % Lower triangular array
        if k >= 2
            % (k , 1:k-1) means entire row k and columns 1 thru k-1
            temp = L(k, 1:k-1);
            L(k, 1:k-1) = L(i, 1:k-1);
            L(i, 1:k-1) = temp;
        end %endif
    end %end outer if ~=
    % get factors L and U
    for j = k+1:n
        % set multiplier in L
        L(j, k) = U(j, k) / U(k, k); % flop count
        % for all in row j = for all in row j - multiplier*for all in row k
        U(j, :) = U(j, :) - L(j, k) * U(k, :)
    end % end gauss elim
end

% solving for c
c = [];
inv(L)
Pb = P*b
c = L\Pb

% solving for x
x = [];
x = U\c
