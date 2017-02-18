%Jacobi Method, Modified:
% NUMERICAL METHODS: MATLAB Programs
%(c) 1999 by John H. Mathews and Kurtis D. Fink
%To accompany the textbook:
%NUMERICAL METHODS Using MATLAB,
%by John H. Mathews and Kurtis D. Fink
%ISBN 0-13-270042-5, (c) 1999
%PRENTICE HALL, INC.
%Upper Saddle River, NJ 07458

%P = zeros(2,1)
%u = matlabFunction((5-v)/3)
%v = matlabFunction((5-u)/2)
%TODO get input
A = [3, 1; 1, 2]
b = [5; 5]
X = zeros(2,1)
iterations = 10;
N = length(b);
Tol = 0.00000001
%TODO, make strictly diagonal

for k=1:iterations
   for j=1:N
      X(j)=(b(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);
   end
   err=abs(norm(X'-P));
   relerr=err/(norm(X)+eps);
   P=X';
   if (err<Tol)|(relerr<Tol)
     break
   end
end