% SOR Successive Over-Relaxation Method
%   This function solves linear equation systems such as Ax=b using SOR 
%   method (Successive Over-Relaxation).
%   When the relaxation scalar w = 1, the method used is Gauss-Seidel.
% Input: A = matrix A 
%        xold = initial x vector 
%        b = solution vector
%        w = relaxation parameter
%        max_it = maximum iterations
%        tol = tolerance

% Output: x = solution
%   error = relative forward error
%   iter = number of iterations
%   boolean flag = 1 means divergence, 0 means convergence
function[x, error, iter, flag]  = SOR_(A, xold, b, w, max_it, tol)

opCount = 0;
% input initialization
% A = [3 1 -1; 
%      2 4 1; 
%     -1 2 5];
% 
% [n, n] = size(A);
% 
% xold = [0; 0; 0];
% b = [4; 1; 1];
% w = 1.25;
% max_it = 2;
% tol = 0.00000001;
x = xold;
diagonalDominant = true;

n = length(A);

if (w == 1)
    fprinf('The relaxation scalar omega = 1. The method used is now Gauss-Seidel')
end

tic;
% check for diagonally dominant convergence guarantee
for r=1:n 
    % Sum of the entire row minus the diagonal
    rowSum = sumabs(A(r,:)) - abs(A(r, r)); 
    opCount = opCount + 1;
    % Check if diagonal is strictly greater than the row sum
    if (abs(A(r,r)) < rowSum)
    % If not, note it
        diagonalDominant = false;           
        break;
    end
end

if (diagonalDominant == false)
% Let user know that convergence is not guaranteed
    fprintf('Matrix A is not strictly diagonal dominant and may not converge.\n');
end

if size(b) ~= size(xold)
    error('The given approximation vector does not match the x vector size');
end


flag = 0;    
iter = 0;
count = 1;
max_it = 2;

% matrix splitting 
% TODO: opCount will be affected by these operations. Need to add.
D = diag(diag(A));
L = tril(A-D);
U = triu(A-D);
%M = D + w*L;
%N = (1 - w)*D - w*U;
%G = M\N;
G = (inv(D+w*L))*(((1-w)*D-w*U)*xold +w*b);
opCount = opCount + 1;
fprintf('Iteration 1: ', G);
datasave = [];
% begin iteration
for iter = 1:max_it                         
        xnew = G;
        RelForError = (norm(xnew-xold))/(norm(xnew));
        opCount = opCount + 1;
        % update approximation
        while (RelForError > tol)
            xold = xnew;
            opCount = opCount + 1;
            G = (inv(D+w*L))*(((1-w)*D-w*U)*xold +w*b)
            xnew = G;
            opCount = opCount + 1;
            RelForError = (norm(xnew-xold))/(norm(xnew));
            if (RelForError <= tol)
                break
            end
            count = count+1;
            x = [x, xnew];
            datasave = [datasave; count, RelForError, flag];
        end
end

    b = b / w % vector b
    if (RelForError > tol) 
       flag = 1;
       fprintf('Did not converge')
   end
        
   % for function return
   x = xnew;
   error = RelForError;
   iter = count;

fprintf('Number of operations: %d\n', opCount);
toc;
fprintf('\n');

fprintf('  iteration    error    flag\n')

disp(datasave)
fprintf(' x final\n')
disp(xnew)
end

