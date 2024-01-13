%This does PCA anallysis of HP data

%clean everything before start
clear all
close all
clc

%%% this is the name of hte hd5file assigned to you
hd5fp = 'pos10.hd5'
hd5fp = '../../../Data/pos10.hd5'

%load the file
h5disp(hd5fp) %diplay hdf5 file to see what's inside
%%inter=h5read(hdfp,'/inter'); %read inter helical parameters
hp=h5read(hd5fp,'/inter'); %read inter helical parameters

%% assign 
% Roll=inter.Roll;
% Tilt=inter.Tilt;
% Shift=inter.Shift;
% Slide=inter.Slide;
% Rise=inter.Rise;
% Twist=inter.Twist;

%% reorganize all the data into a single set that we'll manipulate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CHANGE ONLY THESE LINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
data=vertcat(hp.Roll);
namelist={'Roll'};


%%%%%%%%%%%%%%%%%%%%%%%%%%
datadivs = 4; 
[xdim,tdim] = size(data);
tsub = tdim/datadivs;

subdata = zeros(xdim,tsub);
meandat = zeros(xdim);
covA    = zeros(xdim,xdim);
vec = zeros(xdim,xdim);
val = zeros(xdim,xdim);

figure(1);
pcnt = 0;
for dset = 1:datadivs
    lo=(dset-1)*tsub + 1;
    hi = dset*tsub;
    subdata =data(:,lo:hi);
    meandat = mean(subdata,2);
    covA =cov(subdata');
    [vec,val] =eig(covA);
    %figure(dset);
    pcnt = pcnt + 1; 
    subplot(datadivs,2,pcnt);
    contour(abs(vec'*vec)>0.3);
    pcnt = pcnt + 1; 
    subplot(datadivs,2,pcnt);
    plot(cumsum(flip(diag(val)))*100/sum(diag(val)),'o')
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%HERE's the MATLAB way for the entire exercise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%[coeffs,score,pcs] = pca(data(:,20000:70000)', 'Algorithm','eig');
%% coeffs is similar to vec ... except that vec is reported 
%%   in the opposite order of coeffs. thus coeffs(1) = vecs(end)
%%  pcs is similar to the projection of data onto the vecs.

% figure(1);
% plot(-vec(:,end),coeffs(:,1),".")  %% a plot comparing the PC1 from two methods
% figure(2);
% plot(-vec(:,1),coeffs(:,end),".")  %% a plot comparing PCn from the two methods
% figure(3);
% plot(1:147,-vec(:,1),1:147,coeffs(:,end))  %% a plot of PC1 from the two methods

evals = diag(val)
tmp1 = vec(:,end);
tmp2 = vec(:,end-1);
x = 1:length(evals);

toppercent = cumsum(flip(evals/sum(evals)*100));

figure(4)
plot(x,toppercent,'-O',x,33*ones(1,length(x))) 
title("Evals");
nomodes = sum(toppercent < 33)-1 ;


figure(5)
hold on;
for index = 0:nomodes
    l1 = sprintf('#%02d  %5.1f%%', index+1, evals(end-index)/sum(evals)*100) ;
    plot(x,vec(:,end-index),'DisplayName',l1);
end
title("Evectors");
legend
hold off;

evals = diag(val)
tmp1 = vec(:,end);
tmp2 = vec(:,end-1);
%x = 1:147;

figure(6)
plot(x,tmp1,x,tmp2);
l1 = sprintf(' %8.3f', evals(1)) 
l2 = sprintf('%8.3f', evals(2))
legend(l1,l2)
title("E1 and E");

PC1coeff = (data-meandat)'*tmp1;
pc1data = meandat' + PC1coeff*tmp1';
fakeHP = pc1data';
figure(7)
[X,Y] = meshgrid(-73:1:73,1:21);
for i = 1:1
pltdat = fakeHP((i-1)*147+1:(i)*147,1:21);
%%ax1 =subplot(1,1,plotindex(i));
surf(X,Y,pltdat')
%%heatmap(pltdat')
title(namelist(i))
xlim([-74,74])
end

figure(8)
for i = 1:1
pltdat  = vec((i-1)*147+1:(i)*147,end);
pltdat2 = vec((i-1)*147+1:(i)*147,end-1);
%ax1 =subplot(1,1,plotindex(i));
plot(1:147,pltdat',1:147,pltdat2)
title(ax1,namelist(i))
%xlim([-74,74])
end
