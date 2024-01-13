%% this is the matlab code for assignment 3
%%% Read Toms data
%version 1.0 lla027@latech.edu
%%Read my data

close all, clear all, clc;

hd5fp = strcat('hps.hd5');

h5disp(hd5fp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialize data selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin=15-6; %getting position plus 6 for my assignment.
xmax=xmin+146;
%read inter roll parameters
data=h5read(hd5fp,'/inter').Roll(xmin:xmax,1:1000:100000);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OVERVIEW OF DATA (?) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xdim,tdim]=size(data);%setting dimensions
xvals = 1:xdim; 
meanvals =mean(data,2);
maxvals=max(data,[],2); 
minvals=min(data,[],2);
stdvals=std(data,[],2);
%% make some pretty plots now
%%

data = data - meanvals;
figure(1)
surf(data(:,:))
title("Data - Mean")

%% 
% Covariance matrix & calc of eigval and eig vecs
covA = cov(data'); 
[evecs,evals] =eig(covA,'vector');
% Sorted eigen values
maxi = 147;
maxj = 100;
toppercent = cumsum(flip(evals)/sum(evals));
enumpcnt = 0.25;
enumcnt = sum(toppercent<enumpcnt);

figure(2)
surf(evecs'*evecs)  
%%% columns of evecs are the PCs 
%%   they should be orthonormal 
%% above plot says they are order from small to large (had to flip)
title("orthogonality of evecs")

%% now filter out the coeffs before the recon
%%   set filtercoeff to the coeffs then zero out the ones you don't want
coeffs = data'*evecs;
filtercoeff = coeffs;
filtercoeff(:,1:end-enumcnt) = zeros(maxj,maxi-enumcnt);

figure(3)
%%% remember we have removed the "mean value"
%% so this data is "centered" 
%% have to add meansurf to obtain original data
%% see next slides
subplot(1,2,1)
surf(fliplr(evecs*coeffs'));
title("My Full Recon of Centered")
subplot(1,2,2)
projection = evecs*filtercoeff';
surf(fliplr(projection));
title("Filtered Recon of Centered") 


meansurf = meanvals; 
figure(4)
subplot(1,3,1)
surf(data + meansurf);
title("Raw Data")
subplot(1,3,2)
surf(projection + meansurf)
title("Filtered Recon")
subplot(1,3,3)
surf(projection- data)
title("Residual")
%% 

%%%%%% Filter in 2d %%%%%%%
figure(5)
fftdat2d = fft2(data);
bar(abs(fftdat2d))
%% 
minamp = 1e-4;
ids = abs(fftdat2d) > minamp*xdim*100000;
nfftdat2d = fftdat2d.*ids;
smoothdat = ifft2(nfftdat2d);
figure(6)
subplot(1,2,1)
surf(smoothdat())
title("2D filtered")
subplot(1,2,2)
surf(ifft2(fftdat2d))
title("NO filtered")
%% 

% FFT2D and PCA No filter
figure(7)
subplot(1,2,1)
surf(ifft2(fftdat2d))
title("FFT no filter")
subplot(1,2,2)
surf(fliplr(evecs*coeffs'));
title("My Full Recon of Centered")
%% 

%FFT and PCA after filter
figure(8)
subplot(1,2,1)
surf(smoothdat() + meansurf)
title("2D filtered")
subplot(1,2,2)
surf(projection + meansurf)
title("Filtered Recon")
%% 

figure(9)
subplot(1,2,1)
surf()