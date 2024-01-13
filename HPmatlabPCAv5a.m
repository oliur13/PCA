%This does PCA anallysis of HP data

%%%% CAUTION STILL A WORK IN PROGRESS
%%  NOT SURE THE FINAL CALCS FOR RECONSTRUCTION ARE CORRECT
%%%   tcbishop 2/1/22  USE AT OWN RISK

%clean everything before start
clear all
close all
clc

%%% this is the name of hte hd5file assigned to you
hd5fp = '../../../Data/pos10.hd5'
%%%%

%load the file
h5disp(hd5fp) %diplay hdf5 file to see what's inside
intra=h5read(hd5fp,'/inter'); %read inter helical parameters
data=vertcat(intra.Roll);
namelist={'Roll'};
%% data consists of  
% inter.Roll;
% inter.Tilt;
% inter.Shift;
% inter.Slide;
% inter.Rise;
% inter.Twist;
% intra.Shear
% intra.Stretch
% intra.Stagger
% intra.Buckle
% intra.Propeller
% intra.Opening



%%%%%%%%%%%%%%%%%%%%%%%%%%
datadivs = 4; 
[xdim,tdim] = size(data);
tsub = tdim/datadivs;
%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize variables to start intermediate values
subdata = zeros(xdim,tsub);
meandat = zeros(xdim);
covA    = zeros(xdim,xdim);
vec = zeros(datadivs,xdim,xdim);
val = zeros(datadivs,xdim,xdim);
tmp1 = zeros(xdim,xdim);
nomodes = zeros(datadivs,1);

%% for each subdata set
for dset = 1:datadivs
    lo=(dset-1)*tsub + 1;
    hi = dset*tsub;
    %% select the subdata
    subdata =data(:,lo:hi);
    %% calculate mean and covariance of this subset
    meandat = mean(subdata,2);
    covA =cov(subdata');
    %% compute the eigenvectors/values of covariance
    %% this is the PCA
    [vec(dset,:,:),val(dset,:,:)] =eig(covA);
    
end
%% above is the data analysis loop
%% below is same loop but over figure 1 
%% could have combined them but easier to track this way

figure(1);
title(namelist(1))
 %% Figure 1
 %%  includes vector dot products
pcnt = 0;
hold on;
for dset = 1:datadivs
    %% now make some figures and plots
    pcnt = pcnt + 1; 
    %subplot(datadivs,2,pcnt);
    tmp1 = reshape(vec(dset,:,:),[xdim,xdim]);
   % contour(abs(tmp1'*tmp1));
    pcnt = pcnt + 1; 
    %% and 
    %subplot(datadivs,2,pcnt);
    dvals = diag(reshape(val(dset,:,:),[xdim,xdim]));
    toppercent = cumsum(flip(dvals))*100/sum(dvals);
    plot(toppercent);
    nomodes(dset) = sum(toppercent < 33)-1 ;
end
hold off;


%% make the 4x4 plot of all the cross products
figure(2);
title(namelist(1))
pcnt = 0;
for dset1 = 1:datadivs
    for dset2 = 1:datadivs
      pcnt=pcnt+1;
      tmp1 = reshape(vec(dset1,:,:),[xdim,xdim]);
      tmp2 = reshape(vec(dset2,:,:),[xdim,xdim]);
      subplot(datadivs,datadivs,pcnt);
      %%contour(flipud(abs(tmp1'*tmp2) > 0.5)); 
      h=heatmap(abs(tmp1'*tmp2));
      h.GridVisible = 'off';
      h.ColorLimits = [0 1];
      h.YLabel = sprintf("E%02d",dset1);
      h.XLabel = sprintf("E%02d",dset2);
      ax = gca;
      ax.XDisplayLabels = nan(size(ax.XDisplayData));
      ax.YDisplayLabels = nan(size(ax.YDisplayData));
      if(mod(pcnt,datadivs) )
           h.ColorbarVisible = 'off';
      else
          h.ColorbarVisible = 'on';
      end
       %%% still need to fix the axis numbering y axis goes wrong way
    end
end

figmin=2;
for dset = 1:datadivs
    lo=(dset-1)*tsub + 1;
    hi = dset*tsub;
    %% select the subdata
    subdata =data(:,lo:hi);
    evecs = reshape(vec(dset,:,:),xdim,xdim);
    coeffs = subdata'*evecs; %% 10000x147 coefficients (time by vector)
    figmin=figmin+1;
    figure(figmin);
    surf(coeffs);
    title(sprintf("Coeffs %d",dset));
    
    figmin=figmin+1;
    figure(figmin);
    proj = zeros(xdim,tsub);
    enums = xdim:-1:nomodes(dset);
    for i = enums;
        proj  = proj  + (coeffs(:,i)*evecs(:,i)')';
    end
    %% above is only capturing the fluctuation !?!?
    %% still need to add the mean value of subdata back in
    surf(proj(:,1:10:tsub))
    title(sprintf("Proj %d",dset));
    figmin=figmin+1;
    figure(figmin);
    %plot(1:xdim,mean(subdata,2)-mean(proj,2)) %% 
    surf(subdata-proj) %% subdata -proj should look like mean !?!?!
    title(sprintf("Resid %d",dset));
end



