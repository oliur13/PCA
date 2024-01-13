%% this is the matlab code for assignment 3
%%% Read Toms data
%lla027@latech.edu
%%Read my data

close all, clear all, clc;

hd5fp = strcat('hps.hd5');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialize data selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin=15-6; %getting position plus 6 for my assignment.
xmax=xmin+146;
%read inter roll parameters
data=h5read(hd5fp,'/inter').Roll(xmin:xmax,1:1000:100000);
origdata=data;

[xdim,tdim]=size(data);%setting dimensions
xvals = 1:xdim; 
meanvals =mean(data,2);
%%

data = data - meanvals; %data centering
figure(1)
subplot(1,2,1)
surf(data(:,:))
title("Centered data")
subplot(1,2,2)
surf(origdata)
title("Original Data")

figure(10)
subplot(1,2,1)
surf(data(:,:))
title("Centered data")
view(90,0)
subplot(1,2,2)
surf(origdata)
title("Original Data")
view(90,0)
%% 
%%%%%% PCA %%%%%%%
% Covariance matrix & calc of eigval and eig vecs
covA = cov(data'); 
[evecs,evals] =eig(covA,'vector');
% Sorted eigen values
maxi = 147;
maxj = 100;
toppercent = cumsum(flip(evals)/sum(evals));
enumpcnt = 0.25;
enumpcnt2 = 1;
enumcnt = sum(toppercent<enumpcnt);
enumcnt2 = sum(toppercent<enumpcnt2);

%% now filter out the coeffs before the recon
%%   set filtercoeff to the coeffs then zero out the ones you don't want
coeffs = data'*evecs;
filtercoeff = coeffs;
filtercoeff(:,1:end-enumcnt) = zeros(maxj,maxi-enumcnt);
projection = evecs*filtercoeff';

% Filter out all the PCA components
filtercoeff2 = coeffs;
filtercoeff2(:,1:end-enumcnt2) = zeros(maxj,maxi-enumcnt2);
projection2 = evecs*filtercoeff2';

%%% remember we have removed the "mean value"
%% so this data is "centered" 
%% have to add meansurf to obtain original data

meansurf = meanvals;

figure(2)
subplot(2,3,1)
surf(data);
title("Centered Data")
subplot(2,3,2)
surf(projection + meansurf)
title("Filtered Recon")
subplot(2,3,3)
surf(projection- data)
title("Residual")
subplot(2,3,4)
surf(data+meanvals);
title("Centered Data")
subplot(2,3,5)
surf(projection2+meanvals)
title("All principal components Filtered Recon")
subplot(2,3,6)
surf(data-projection2)
title("Residual")
%% 

%%%%%%FFT Analysis in 2d %%%%%%%
% Compute the FFT coefficients of the signal
fftdat = fft2(data);
% Compute the number of data points
N = length(data);
% Initialize a matrix to store the basis functions
basis = zeros(N);
% Extract the basis functions
for k = 1:N
    % Create a vector with a single non-zero element
    x = zeros(N, 1);
    x(k) = 1;
    % Compute the FFT of the vector to obtain the basis function
    basis(:, k) = fft(x);
    % Normalize the basis functions
    % Compute the norm of the k-th column
    norm_k = sqrt(sum(abs(basis(:, k)).^2));
    % Normalize the k-th column by dividing by its norm
    basis(:, k) = basis(:, k) / norm_k;
end

% Compute the matrix of dot products (i.e., the inner products)
orthmat = basis' * basis;
% Display the result as a surface plot
% figure(5)
% surf(abs(orthmat))
%% 
 
minamp = 1e-4;
ids = abs(fftdat) > minamp*xdim*100000;
nfftdat2d = fftdat.*ids;
smoothdat = ifft2(nfftdat2d);

% Filter out all the FFT components
minamp2 = 0;
ids2 = abs(fftdat) > minamp2*xdim*100000;
nfftdat2d2 = fftdat.*ids2;
smoothdat2 = ifft2(nfftdat2d2);

figure(3)
subplot(2,3,1)
surf(data)
title("Raw data")
subplot(2,3,2)
surf(smoothdat)
title("Filtered data")
subplot(2,3,3)
surf(data-smoothdat)
title("Residual")
subplot(2,3,4)
surf(data+meanvals)
title("Raw data")
subplot(2,3,5)
surf(smoothdat2+meanvals)
title("All FFT components filtered data")
subplot(2,3,6)
surf(data-smoothdat2)
title("Residual")
%% 

% FFT2D and PCA all components filtered
figure(4)
subplot(2,3,1)
surf(data+meanvals)
title("Centered data")
subplot(2,3,2)
surf(smoothdat2+meanvals)
title("All FFT components filtered data")
subplot(2,3,3)
surf(data-smoothdat2)
title("Residual")
subplot(2,3,4)
surf(data+meanvals);
title("Centered Data")
subplot(2,3,5)
surf(projection2+meanvals)
title("All principal components Filtered Recon")
subplot(2,3,6)
surf(data-projection2)
title("Residual")

figure(12)
subplot(2,2,1)
surf(data+meanvals)
title("Centered data")
view(90,0)
subplot(2,2,2)
surf(smoothdat2+meanvals)
title("All FFT components filtered data")
view(90,0)
subplot(2,2,3)
surf(data+meanvals);
title("Centered Data")
view(90,0)
subplot(2,2,4)
surf(projection2+meanvals)
title("All principal components Filtered Recon")
view(90,0)
%% 

%FFT and PCA after filter
figure(5)
subplot(1,2,1)
surf(smoothdat + meansurf)
title("FFT filtered Recon with the highest amplitudes")
subplot(1,2,2)
surf(projection + meansurf)
title("PCA filtered Recon with the largest evals")

figure(11)
subplot(2,2,1)
surf(smoothdat + meansurf)
title("FFT filtered Recon with the highest amplitudes")
view(0,0)
subplot(2,2,2)
surf(projection + meansurf)
title("PCA filtered Recon with the largest evals")
view(0,0)
subplot(2,2,3)
surf(smoothdat + meansurf)
title("FFT filtered Recon with the highest amplitudes")
view(90,0)
subplot(2,2,4)
surf(projection + meansurf)
title("PCA filtered Recon with the largest evals")
view(90,0)
%% 
% Orthogonality

figure(6)
subplot(1,2,1)
surf(abs(orthmat))
title('orthogonality of FFT components')
subplot(1,2,2)
surf(evecs'*evecs)
title('orthogonality of evecs')
%% 
% Orthogonality between FFT components and PCA components

figure(7)
surf(abs(basis)*evecs)
title('orthogonality between FFT components and PCA evecs')