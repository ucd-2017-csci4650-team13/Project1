%Jacobi Method, Modified:

% NUMERICAL METHODS: MATLAB Programs
%(c) 1999 by John H. Mathews and Kurtis D. Fink
%To accompany the textbook:
%NUMERICAL METHODS Using MATLAB,
%by John H. Mathews and Kurtis D. Fink
%ISBN 0-13-270042-5, (c) 1999
%PRENTICE HALL, INC.
%Upper Saddle River, NJ 07458

% Function takes in
% NxN Coefficient Matrix A
% Nx1 Solution Matrix b
% Nx1 Matrix of inital guess P
% Returns a matrix of approximate solutions X

function X = Jacobi(augA, P)
rows = size(augA,1)
columns = rows + 1
A = augA
A(:,columns) = []
%A = input('enter an initial matrix: \n');
b = augA(:,columns)
opCount = 0

X = zeros();
iterations = 10;
N = length(b);
Tol = 0.00000001;
diagonalDominant = true;

% Check if A is strictly diagonally dominant
% For each row in matrix A
for r=1:N
    rowSum = sumabs(A(r,:)) - abs(A(r, r)) % Sum of the entire row minus the diagonal
    % Check if diagonal is strictly greater than the row sum
    if abs(A(r,r)) < rowSum
        diagonalDominant = false;           % If not, note ite
        break;
    end
end

if diagonalDominant == false
    fprintf('Matrix A is not strictly diagonal dominant and may not converge.\n');
end

for k=1:iterations
    for j=1:N
        X(j)=(b(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);     % slicing using coefficients to solve for uk and vk
    end
    err=abs(norm(X'-P));
    relerr=err/(norm(X)+eps);
    P=X';
    if (err<Tol)||(relerr<Tol)
        break
    end
end
X = X'

