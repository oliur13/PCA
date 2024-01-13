%%% some tests of qr and other routines

%% we'll use the lehmer matrix from matlab
A = gallery('lehmer',10);

%% call to myqr to get eigenvalues using our GS routine
%% and a orthogonal basis for the matrix
%% this does not compute eigenvectors just values
[q,r,B] = myqr(A,400);

%% check that the resutls are really orthogonal
disp("Verify that q is orto-normal")
tol = 1e-12 ;
%%abs(q'*q - eye(10,10)) < tol 
%% above should be all true
orthotest = q'*q

%%% now run the QR_Decomposition using 
% %%% matlabs qr routine
%%% 
[evals,notEVects ] = QR_Decomposition(A);

%% now get the eigen vectors and eigenvalues 
%%%from matlabs eig
[vecs, vals] = eig(A);

vals = flip(diag(vals)); %%% because in diff orders

%% show that the vals from matlabs eig 
%% and the evals from QR_Decomp are the same
vals
evals

%% use the power method to get 1st evec from power
[a,b] = power_method(A,rand(10,1));

%% here is the largest eigenvector by power_method and by matlabs eig
figure(1)
plot(vecs(:,end),b,"x")
title("Evec and PowerMethod")
figure(2)
plot(notEVects(:,1),b,"o") 
title("Evec and QR")


%% should be a line with slope of 1 
% %%%if they are the sam




