%%% data -> cov -> eig -> PCs to
%%% matlabs PCA routine for doing the same
%%% apparently there are a bunch of switches to pay attention to
%%% that are required to get the same answer for both
clear;
%fpn = ["posp10" ]
%basedir='../../../../rnahome/tmbshare/public_html/sims/601/WT/'
hd5fp = 'hps.hd5'
h5disp(hd5fp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialize data selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin=15+5 %getting position plus 2 for my assignment.SHOULD THIS BE 177?
xmax=xmin+146;
%read inter helical parameters
data=h5read(hd5fp,'/inter').Roll(xmin:xmax,1:10000);
meanvals =mean(data,2);
meansurf =data - meanvals;
%% CENTER THE DATA for all values x
%% centering the data appears to be one of the first issues
%% cov by definition should center the data
%% but PCA and cov seem to be using diff methods to achieve this
%% so remove the mean value (along "x") for all observations (along"t")
data = data - meanvals;
%% data is now centered... which means we'll have to add the meansurf
%% when do the reconstruction
surf(data)
title("Centered Data")
%% THE MATLAB PCA routine
%[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(X,eig)
% [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED,mu]=pca(data','Algorithm','eig','Centered','off','Economy','off');
%% tell PCA to use EIG for eigenvector/value calc. to NOT center the data and to
%% calculate ALL eigenvectors... (economy off)
%%% following is taken from matlab manual
%% rows of matrix given to PCA should be the observations so data(x,t)' = data(t,x)
%% columns of data' are the observations as read from the HD5 file
%% Each column of COEFF contains coefficients for one principal
%% component. The columns are in descending order in terms of component
%% variance (LATENT)
%%% score is representation of X in principal component space
%%% the centered
%%%%%%%%%%%%%%%%%%%
%%%% input data is reconstructed by SCORE*COEFF'
%% data' = SCORE*COEFF' so it follows that
%%%% data' * COEFF = SCORE b/c COEFF' * COEFF = I
%%%%%%%%%%%%%%%%%%%%%%
%% LATENT are the eigenvalues of covariance matrix
%%% COEFF is maxi x maxi
%%% SCORE is maxj x maxi
%%% LATENT is maxi x 1
%%% Explained is maxi x 1
%%% DOING ALL OF THIS MANUALLY
%%%
covA = cov(data'); %% RTFM this should remove the mean value from ea. column
[evecs,evals] =eig(covA,'vector'); %% put evals in a vector instead of matrix
%% columns of evecs are the principal components same as COEFF but flipped LR
%% COEFF(1:10,1) = evecs(1:10,end) .. last evec corresponds largest eval
%% evals are smallest to largest, LATENT is largest to smallest
%% so LATENT and flip(evals) are the same
%% demonstrate that evecs and COEFF are same by taking dot product
figure(1)
surf(abs(evecs'*fliplr(COEFF))) %%% WHATS THIS ALL ABOUT (see line 73)
title("Eig Vecs vs PCA vecs")
%% note that flipped one of them left right b/c they are
%% in different order
figure(2)
plot(mean(data(:,:),2));
title("mean of centered data")
figure(3)
plot(evals)
title("evals from Eig")
%% note they go from smallest to largest
maxi = 147;
maxj=10000;
figure(4)
toppercent = cumsum(flip(evals)/sum(evals));
enumpcnt = 0.28;
enumcnt = sum(toppercent<enumpcnt);
plot(1:maxi,toppercent,'o',1:maxi,enumpcnt*ones(maxi,1))
title("Eval contrib.")
% figure(4)
%
% %%% columns of evecs are the PCs
% %% they should be orthonormal
% %% above plot says they are order from small to large (had to flip)
% title("orthogonality of evecs")
% figure(5)
% plot(1:maxi,evecs(:,end),1:maxi,evecs(:,end-1))
figure(7)
%% coeffs here are same as PCA's score.. flipped LR
%% coeffs here is not PCAs COEFF
coeffs = data'*evecs; %% maxj x maxi (time x pca)
subplot(1,3,1)
surf(coeffs)
title("Coeffs")
subplot(1,3,2)
surf(fliplr(SCORE))
title("PCA Score");
%% now filter out the coeffs before the recon
%% set filtercoeff to the coeffs then zero out the ones you don't want
filtercoeff = coeffs;
filtercoeff(:,1:end-enumcnt) = zeros(maxj,maxi-enumcnt);
subplot(1,3,3)
surf(filtercoeff);
title("Filtered Coeffs")
projection = evecs*filtercoeff';
figure(11)
% subplot(1,3,1)
% surf(meansurf);
% title("Raw Data")
subplot(1,2,1)
surf(projection)
title("Filtered Data")
subplot(1,2,2)
surf(projection - meansurf)
title("PCA Residual")
%%%%%% Filter in 2d %%%%%%%
[xdim tdim]=size(data)%setting dimensions
xvals = 1:147; %%% FIX THIS WHWHGHGHGH
meanvals =mean(data,2);
maxvals=max(data,[],2); %M = max(A, [], 'all'); %taken from reference section
minvals=min(data,[],2);
stdvals=std(data,[],2);
shiftdata = data - meanvals;
figure(5)
fftdat2d = fft2(shiftdata);
fftdat=fft(shiftdata);
minamp = 0.000548787
ids = abs(fftdat2d) > minamp*xdim*100000;
nfftdat = fftdat2d.*ids;
smoothdat = ifft2(nfftdat);
%surf(log10(abs( smoothdat(:,1:500:100000) - shiftdata(:,1:500:100000))))
% subplot(1,3,1)
% surf(meansurf)
% title('Centered Data')
subplot(1,2,1)
surf(smoothdat)
title("Filtered Data")
subplot(1,2,2)
surf(smoothdat-meansurf)
title("FFT Residual")
figure(21)
surf(abs(fftdat2d)/xdim)
title("fourier spectra")
figure(22)
fids = 1:(xdim+1)/2;
plot(fids,ones(length(fids),1)*minamp,fids,sort(abs(fftdat2d(fids))/xdim),'-x')
title('sorted normalize spectra')
figure(23)
subplot(1,2,1)
surf(smoothdat)
title("Filtered FFT")
subplot(1,2,2)
surf(projection)
title("Filtered PCA")
% figure(21)
% bar(abs(fftdat)/xdim)
% title("fourier spectra")
%
% figure(22)
% fids = 1:(xdim+1)/2;
% plot(fids,ones(length(fids),1)*minamp,fids,sort(abs(fftdat(fids))/xdim),'-x')
% title('sorted normalize spectra')
% times vector
N = 100;
dt = 0.5;
tmin = 0.0;
t = tmin + dt*[0:N-1]';
tmax = tmin + dt*(N-1);
% frequency spacing
df = 1/(N*dt);
% number of unknowns same as data
M = N;
% set up least-squares G matrix
G = zeros(N,M);
G(:,1)=ones(N,1);
for p = 2*[1:M/2-1]
G(:,p) = cos(pi*p*df*t);
G(:,p+1) = sin(pi*p*df*t);
end
p=M/2;
G(:,M) = cos(2*pi*p*df*t);
G_norm=bsxfun(@rdivide,G,sqrt(sum(G.^2,1)));
figure(6)
subplot(1,2,1)
surf(G_norm'*G_norm)
title("Orthogonality of FFT")
subplot(1,2,2)
surf(evecs'*evecs)
title("Orthogonality of PCA")