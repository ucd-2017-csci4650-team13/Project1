function [X, errFlag] = Jacobi(augA, P, Tol, iterations, handles)
rows = size(augA,1);
columns = rows + 1;
A = augA;
A(:,columns) = [];
b = augA(:,columns);
opiCount = 0;

X = zeros();
N = length(b);
errFlag = false;
convergenceStr = '';
% Check if A is strictly diagonally dominant
% For each row in matrix A
for r=1:N
    rowSum = sumabs(A(r,:)) - abs(A(r, r)); % Sum of the entire row minus the diagonal
    % Check if diagonal is strictly greater than the row sum
    if abs(A(r,r)) < rowSum
        convergenceStr = 'Matrix A is not strictly diagonal dominant and may not converge.';
        set(handles.lSysOutputText, 'string', convergenceStr);
        break;
    end
end

tic;

for k=1:iterations
    for j=1:N
        X(j)=(b(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);     % slicing using coefficients to solve for uk and vk
    end
    opiCount = opiCount + N^2;
    err=abs(norm(X'-P));
    relerr=err/(norm(X)+eps);
    P=X';
    if (err<Tol)||(relerr<Tol)
        break
    end
end

if errFlag == false
    time = toc;
    opsString = ['Number of Iterations = ', num2str(k), ',Operations = ', num2str(opiCount)];
    timeString = ['Elapsed Time = ', num2str(time), ' seconds'];
    newString = combineString(convergenceStr, opsString);
    newString = combineString(newString, timeString);
    set(handles.lSysOutputText, 'string', newString);
end