function[x, errFlag]  = SOR(augA, x, w, max_it, tol, handles)
errFlag = false;
rows = size(augA,1);
columns = rows + 1;
A = augA;
A(:,columns) = [];
b = augA(:,columns);

iCount = 0;
GSStr = '';
convStr = '';

x = x';

n = length(A);

if (w == 1)
    GSStr = 'The relaxation scalar omega = 1. The method used is now Gauss-Seidel';
    set(handles.lSysOutputText, 'string', GSStr);
end

tic;
% check for diagonally dominant convergence guarantee
for r=1:n
    % Sum of the entire row minus the diagonal
    rowSum = sumabs(A(r,:)) - abs(A(r, r));
    iCount = iCount + 1;
    % Check if diagonal is strictly greater than the row sum
    if (abs(A(r,r)) < rowSum)
        % If not, note it
        fprintf('Matrix A is not strictly diagonal dominant and may not converge.\n');
        convStr = 'Matrix A is not strictly diagonal dominant and may not converge';
        newStr = combineString(GSStr, convStr);
        set(handles.lSysOutputText, 'string', newStr);
        break;
    end
end

if size(b) ~= size(x)
    errStr = 'The given approximation vector does not match the x vector size';
    set(handles.lSysOutputText, 'string', errStr);
    errFlag = true;
else
    flag = 0;
    iCount = 1;
    % matrix splitting
    % TODO: opiCount will be affected by these operations. Need to add.
    D = diag(diag(A));
    L = tril(A-D);
    U = triu(A-D);
    %M = D + w*L;
    %N = (1 - w)*D - w*U;
    %G = M\N;
    G = (inv(D+w*L))*(((1-w)*D-w*U)*x +w*b);
    iCount = iCount + 1;
    %fprintf('Iteration 1: ', G);
    datasave = [];
    % begin iteration
    for iter = 1:max_it
        xnew = G;
        RelForError = (norm(xnew-x))/(norm(xnew));
        iCount = iCount + 1;
        % update approximation
        while (RelForError > tol)
            x = xnew;
            iCount = iCount + 1;
            G = (inv(D+w*L))*(((1-w)*D-w*U)*x +w*b);
            xnew = G;
            iCount = iCount + 1;
            RelForError = (norm(xnew-x))/(norm(xnew));
            if (RelForError <= tol)
                break
            end
            x = [x, xnew];
            datasave = [datasave; iCount, RelForError, flag];
        end
    end
    
    b = b / w; % vector b
    if (RelForError > tol)
        flag = 1;
        convStr = 'Did not converge';
        set(handles.lSysOutputText, 'string', convStr);
    end
    
    % for function return
    x = xnew;
    error = RelForError;
    iter = iCount;
    
    time = toc;
    opsStr = ['Number of Iterations = ', num2str(iter), ',Operations = ', num2str(iCount)];
    timeString = ['Elapsed Time = ', num2str(time), ' seconds'];
    newStr = combineString(GSStr, convStr);
    newStr = combineString(newStr, opsStr);
    newStr = combineString(newStr, timeString);
    set(handles.lSysOutputText, 'string', newStr);
    x=x';
    %     fprintf('  iteration    error    flag\n')
    %
    %     disp(datasave)
    %     fprintf(' x final\n')
    %     disp(xnew)
end