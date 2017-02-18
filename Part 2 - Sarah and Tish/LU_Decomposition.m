clc;
clear all;
close all;

% matrix A
A=[-1, 0, 1;
   2, 1, 1;
   -1, 2, 0]

% size of n x n
[n, n] = size(A)
col = n

% init U to A
U = A

% init L to identity
L = eye(n)

% identity matrix for pivoting
P = eye(n)

% row contining max value
max = 0;
row = 0;
for i = 1:n
    if abs(max) < abs(A(i,1))
        max = A(i,1)
        row = i
    end
end

if row > 1
    %interchange row with 1
    temp = A
    temp_p = P
    for i = 1:n
        % move row 1 into the row with max pivot 
        A(row, i) = A(1, i)
        % permute Identity matrix
        P(row, i) = P(1, i)
        % move max pivot row into row 1 
        A(1, i) = temp(row, i)
        % permute Identity matrix
        P(1, i) = temp_p(row, i)
    end
end

U = A  
% Upper triangular matrix
for j = 1:n-1 % n-1 = num of rows - the 1st row
    % check for division by zero
    if U(j,j) == 0
        error('zero pivot encountered');
    end
    %%%%%%%%%%%%%%DEBUG%%%%%%%%%%%%%%%%%
    % eliminate col j to put zero in each location below diag
    % ex. a(j+1, j), a(j+2, j),...a(n,j)
    % var to keep track of val of pivot
    row = j
    for i = j+1:n
        if abs(U(i, 1)) > abs(U(row, 1))
        %interchange row i with row j
        temp_u = U
        temp_p = P
            for index = 1:n
            % move row j into the row i 
            U(row, index) = U(i, index);   
            % permute Identity matrix
            P(row, index) = P(i, index);
            % move row i into row j 
            U(i, index) = temp_u(row, index)
            % permute Identity matrix
            P(i, index) = temp_p(row, index)
            end
        end
        % row multiplier
        % TODO A(j,j) == 0 fix
        multi = U(i, j) / U(j, j)
        
        % L factor gets multiplier
        L(i, j) = multi   
        
        % subtract multiplier * the row from 
        for index = j:col
            U(i, index) = U(i, index) - (multi * U(j, index))
            U(i, index)
        end
        %%%%%%%%%%%%%%BREAKPOINT%%%%%%%%%%%%%%%%%
        row = row + 1
    end
end





