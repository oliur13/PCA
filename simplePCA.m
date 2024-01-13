clear;
maxi = 147;
maxj = 1000;

data = zeros(maxi,maxj);
for i = 1:maxi
for j = 1:maxj
data(i,j) = sin(2*pi*(mod(i,14)*i/maxi))*cos(2*pi*i*j/maxj) + cos(2*pi*8*i/maxi);
end
end

meansurf = repmat(mean(data,2),1,maxj);
%% CENTER THE DATA for all values x
data = data - meansurf;
surf(data)
title("Centered Data")
%% THE MATLAB  PCA routine
%[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(X,eig)
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(data','Algorithm','eig','Centered','off','Economy','off');


%% rows of X are observations
%% columsn of X are variables so X is data'
%% Each column of COEFF contains coefficients for one principal
%%    component. The columns are in descending order in terms of component
%%    variance (LATENT)
%%% score is representation of X in principal component space
%%%      the centered data is reconstructed by SCORE*COEFF'
%% LATENT are the eigenvalues of covariance matrix
%%% COEFF is maxi x maxi 
%%% SCORE is maxj x maxi
%%% LATENT is maxi x 1
%%% Explained is maxi x 1

%%% DOING ALL OF THIS MANUALLY
%%%
covA = cov(data'); %% RTFM this removes the mean value in both dimensions
[evecs,evals] =eig(covA,'vector'); %% put evals in a vector
%% columns of evecs are the principal components same as COEFF but flipped LR
%%      COEFF(1:10,1) = evecs(1:10,end) .. last evec corresponds largest eval
%% evals are smallest to largest, LATENT is largest to smallest
%%    so  LATENT and flip(evals) are the same

%% demonstrate that evecs and COEFF are same 
%%surf(abs(evecs'*fliplr(COEFF)))

figure(2)
plot(mean(data,2));
title("mean")

figure(3)
plot(evals)
title("evals")
figure(4)
toppercent = cumsum(flip(evals)/sum(evals));
plot(toppercent,'o')
enumcnt = sum(toppercent<.33);
title("eval contrib")

figure(6)
contour(evecs'*evecs)  
%%% columns of evecs are the PCs 
%%   they should be orthonormal 
%% above plot says they are order from small to large (had to flip)
title("orthoc of evecs")

% figure(5)
% plot(1:maxi,evecs(:,end),1:maxi,evecs(:,end-1))

figure(7)
%% coeffs here are same as PCA's score.. they are not COEFF
coeffs = data'*evecs; %% 1000x147 (time x pca)
subplot(1,3,1)
surf(coeffs)
title("Coeffs")
subplot(1,3,2)
surf(fliplr(SCORE))
title("PCA Score");
%% now filter out the coeffs before the recon
filtercoeff = coeffs;
%%enumcnt = 50;
filtercoeff(:,1:end-enumcnt) = zeros(maxj,maxi-enumcnt);
subplot(1,3,3)
surf(filtercoeff);
title("Filtered Coeffs")


figure(10)
subplot(1,3,1)
surf(fliplr(evecs*coeffs'));
title("My Full Recon of Centered")
subplot(1,3,2)
surf(SCORE*COEFF')
title("Matlab Full Recon of Centered")
subplot(1,3,3)
projection = evecs*filtercoeff';
surf(fliplr(projection));
title("Filtered Recon of Centered") 


figure(11)
subplot(1,3,1)
surf(data + meansurf);
title("Raw Data")
subplot(1,3,2)
surf(projection + meansurf)
title("Filtered Recon")
subplot(1,3,3)
surf(projection- data)
title("Residual")


figure(12)
subplot(1,2,1)
plot(mean(data,2))
title("Mean of Raw Data")
subplot(1,2,2)
plot(mean(projection,2))
title("Mean of Filtered")

