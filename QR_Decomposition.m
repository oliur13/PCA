%determines eigenvalues using QR Decomposition
%with Q and R calcualted by Matlab
%% 

function [e,A]=QR_Decomposition(A)
e=1;
tol=1e-10;
while(1)
    e_old = e;
    [Q,R]= qr(A); %% using matlabs qr decomposition
    A = R*Q;
    e = diag(A);
     if abs(max(e)-max(e_old)) < tol
         break;
     end
end
end