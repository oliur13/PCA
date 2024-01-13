%% this is the matlab code for assignment 3
%%% Read Toms data
%%%Group Matlab code for assignment3
%version 1.0 lla027@latech.edu
%%Read my data


clear all
close all
clc


hd5fp = 'hps.hd5'

h5disp(hd5fp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialize stuff (???? ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read inter helical parameters
inter=h5read(hd5fp,'/inter');
% reorganize all the data into a single set that we'll manipulate
data=vertcat(inter.Roll);
figure(751)
plot(data(:,500))
%%setting up our specific position

xmin=15-6 %getting position plus 2 for my assignment.SHOULD THIS BE 177?
xmax=xmin+146
data=data(xmin:xmax,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OVERVIEW OF DATA (?) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(231)
plot(data(:,500))
title("t=500")
namelist={ 'Roll'}
%%% calc avg
%%Average values
[xdim ydim]=size(data)%setting dimensions
xvals=1:xdim; %changef to 176
xrange = xvals;
norm = range(xrange);
meanvals =mean(data,2);
maxvals=max(data,[],2); %M = max(A, [], 'all'); %taken from reference section
minvals=min(data,[],2);
stdvals=std(data,[],2);
%% make some pretty plots now
%%
figure(7)
plot(meanvals)
title("mean Roll vals")
% plot of max, min, mean, std deviation
figure(1)
plot(xvals,maxvals,xvals,minvals,xvals,meanvals,xvals,meanvals+stdvals, ...
   xvals,meanvals-stdvals)
   axis tight

   title("Roll Stats",'FontSize',14)
   xlabel('x','FontSize',12)
   ylabel('Value','FontSize',12)
   legend("Max","Min","Mean","+1 \sigma","-1 \sigma",'Location','best')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manual Fitting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the function F(x)
n = 14;
xrnge = xmax - xmin;
x = linspace(-xrnge/2, xrnge/2, 147);
x=x';
F = @(x,C) C(1) + C(2)*x + C(3)*x.^2 + C(4)*sin(2*pi*n*x/xrnge) + C(5)*cos(2*pi*n*x/xrnge);

% Define the objective function to be minimized
obj = @(C) sqrt(sum((meanvals - F(x,C)).^2))/xdim;

% Find the optimal values for the constants using fminsearch
C0 = [0 0 0 0 0]; % initial guess for the constants
C = fminsearch(obj,C0);

% Report the optimal values for the constants
disp(['C1 = ' num2str(C(1))])
disp(['C2 = ' num2str(C(2))])
disp(['C3 = ' num2str(C(3))])
disp(['Cs = ' num2str(C(4))])
disp(['Cc = ' num2str(C(5))])

% Plot the mean value of the Roll data and F(x)
figure(999)
plot(x, meanvals, '-b')
hold on
plot(x, F(x,C), '-r')
xlabel('x')
ylabel('Value')
title('Manual fitting for mean values of Roll data')
legend('Mean value of Roll data', 'F(x)')
%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fun with FFTs (FFT analysis of data?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute FT of the avg using  matlabs FFT 
%%fast fourier transform
fftdat=fft(meanvals);
minamp = 0.75;

%%%plot(abs(fftdat(1:(xdim+1)/2))/xdim)
%%%  plots spectra ??? RTFM (read the free manual to figure out how MATLAB does this)
%fft
figure(21)
bar(abs(fftdat))
title("fourier spec")
figure(22)
fids = 1:(xdim+1)/2;
plot(fids,ones(length(fids),1)*minamp,fids,sort(abs(fftdat(fids))/xdim),'-x')
title('sorted normalize spectra')
figure(23)
plot((abs(fftdat(1:(xdim+1)/2))/xdim),'-x')
title('normalized spectra')
%%% pick wave numbers based on spectra  (HINT:  [a,b] = sort (something)  … again RTFM)
[A,B]=sort(abs(fftdat(1:(xdim+1)/2))/xdim)
%%% compute invFFT and show that the FILTERED data is within some specified range of the original data
figure(24)
bar(abs(fftdat.*(fftdat > minamp*xdim))/xdim)
title("filtered spectra")
%% 

ids = abs(fftdat) > minamp*xdim;
nfftdat = fftdat.*ids;

figure(100)
%plot(xvals,real(ifft(nfftdat)))

plot(xvals,meanvals,...
    xvals,real(ifft(nfftdat(:,1))),...
    xvals,minvals,':b',...
    xvals,maxvals,':b')
legend("Original Avarage Values","Inverse Fourier Transform","Minvals","Maxvals")

%%%%%% now in 2d %%%%%%%
% fftdat2d = fft2(data);
% minamp = 0.1
% ids = abs(fftdat2d) > minamp*xdim*100000;
% nfftdat2d = fftdat2d.*ids;
% smoothdat = ifft2(nfftdat2d);
% figure(998)
% surf(smoothdat(:,1:100:100000))