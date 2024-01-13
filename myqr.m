%%% calculte eigenvalues by qr decomposition
%%% does not calculate eigenvectors
%% uses a Gram-Schmidt process to  obtain q
%% Bishop 1/2022
function [q,r,A] = myqr(A,ncnt)
tol = 1e-6;
% ncnt = 0;
maxcnt = 1000;
%% initialize 
[n m]=size(A);  % n= # of rows in A, m= # of columns in A
if n ~= m    % check if A is a square matrix
    d ='ERROR' ;   % if not, display error message
    disp('The matrix has to be square');
else
    ndim = n;
end
q = zeros(n,n);
r = zeros(n,n);

%%%% Gram-Schmidt process
for icol = 1:ndim;
%        fprintf(".") ;
       q(:,icol) = A(:,icol);
%        if(icol==1) 
%            disp(q(:,icol));
%        end
       for kcol = 1:icol-1;
%            fprintf("x");
           r(kcol,icol) = r(kcol,icol) - dot(q(:,kcol),A(:,icol));
           q(:,icol) = q(:,icol) - dot(q(:,kcol),A(:,icol))*q(:,kcol);
%            disp(r);
       end
       r(icol,icol) = norm(q(:,icol),2);
       q(:,icol) = q(:,icol)/norm(q(:,icol),2);
end
%%%  a debugging/validation step that's not necessary
%%fprintf("\n");
%% can now test if q*q' = eye(n,n)
% if ( abs(q*q'- eye(ndim,ndim)) < 1e-14 ) 
%     disp(" Q passed")
% else
%     disp(q*q');
% end

%%% now start all the interative multiplication of rq
%%% to find the eigenvalues (not eigenvectors)
    A = r*q;
    if( all(abs(tril(A,-1) < tol)) ) %% if converged to upper triang
        disp("converged");
        disp(ncnt);
        disp(A);
        return;
    else
        ncnt = ncnt + 1; %% if not exceeded max count
        if(ncnt > maxcnt ) 
            disp(ncnt)
            disp(A)
            return;
        else
%             fprintf(" CNT: %d \n", ncnt);
%             disp(A)
            myqr(A,ncnt);
        end
    end
return  

