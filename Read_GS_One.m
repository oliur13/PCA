%This code reads the hdf5 file, parses the data,
% and makes some plots
% TcB 1/21/22
%clean everything before start
clear all
close all
clc

%%% this is the name of hte hd5file assigned to you
hd5fp = './pos10.hd5'
%%%%

%load the file
h5disp(hd5fp) %diplay hdf5 file to see what's inside
inter=h5read(hd5fp,'/inter'); %read inter helical parameters
%intra=h5read('pos11.hd5','/intra'); %read intra helical parameters
%%   there's more here than we'll use to ignore the rest.

%% assign 
% Roll=inter.Roll;
% Tilt=inter.Tilt;
% Shift=inter.Shift;
% Slide=inter.Slide;
% Rise=inter.Rise;
% Twist=inter.Twist;

%% Pick the single data set that we'll manipulate
data=vertcat(inter.Roll);
namelist={'Roll'};

%% set up some info for plotting etc.

xvals=-73:1:73;
avgvals = zeros(147,1);  %% will hold mean value for each of the 6 sets
stdval  = zeros(1);
maxval  = zeros(1);
minval  = zeros(1);

mintime = 1;
maxtime = 100000;%% use all data
    
    data2=data(:,mintime:maxtime);
    avgvals=reshape(mean(data2,2),147,1);
    datacov = cov(data2);
    %%% we'll extract 

for i = 1:1
   
    nostds = 2;
   
    stdval(i) = std(avgvals(:,i));
    maxval(i) = max(avgvals(:,i));
    minval(i) = min(avgvals(:,i));
    
    stdP = (  nostds*stdval(i) + mean(avgvals(:,i)))*ones(147) ;
    stdM = ( -nostds*stdval(i) + mean(avgvals(:,i)))*ones(147) ; 

    figure(1)
    subplot(1,1,i)
    plot(xvals,avgvals(:,i),'-k',xvals,stdP,'--r', xvals, stdM,'--r');
    
    figure(2)
    subplot(1,1,i)
    surf(data((i-1)*147+1:i*(147),mintime:1000:maxtime))
    axis tight
    plotname=namelist(i);
    title(plotname,'FontSize',14)
    xlabel('x','FontSize',12)
    ylabel('Value','FontSize',12)
end

