%% this is the matlab code for assignment 3
%%% Read Toms data
%%%Group Matlab code for assignment3
%version 1.0 tkm020@latech.edu
%preparing for 2/23/23 lecture
%%Read my data
%fpn = ["hps.hd5" ]
%basedir='C:\Files\PHYS 540\02_06'
hd5fp = strcat('hps.hd5')

h5disp(hd5fp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialize data selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin=15-6 %getting position plus 2 for my assignment.SHOULD THIS BE 177?
xmax=xmin+146;
%read inter helical parameters
data=h5read(hd5fp,'/inter').Roll(xmin:xmax,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OVERVIEW OF DATA (?) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xdim tdim]=size(data)%setting dimensions
xvals = 1:xdim; %%% FIX THIS WHWHGHGHGH
meanvals =mean(data,2);
maxvals=max(data,[],2); %M = max(A, [], 'all'); %taken from reference section
minvals=min(data,[],2);
stdvals=std(data,[],2);
%% make some pretty plots now
%%
figure(2)
plot(xvals,meanvals)
title("mean Roll vals")
% plot of max, min, mean, std deviation
figure(3)
plot(xvals,maxvals,':r',...
    xvals,minvals,':r',...
    xvals,meanvals,'-b')
   axis tight
   title("Roll Stats",'FontSize',14)
   xlabel('x','FontSize',12)
   ylabel('Value','FontSize',12)
   legend("Max","Min","Mean")

figure(1)
surf(data(:,1:500:end))
title("Data(:,1:500:end)")

shiftdata = data - meanvals;
figure(4)
surf(shiftdata(:,1:500:end))
title("Data(:,1:500:end) - Mean")


%%%%%% Filter in 2d %%%%%%%
figure(5)
fftdat2d = fft2(shiftdata);
minamp = 0.1
ids = abs(fftdat2d) > minamp*xdim*100000;
nfftdat2d = fftdat2d.*ids;
smoothdat = ifft2(nfftdat2d);
surf(smoothdat(:,1:500:100000))
title("2D filtered")