%ShowPosterior
%   Plot posterior distributions from MCMC samples
%
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

clear
close all
clc

Wi = 360;
Le = 330;

ftsz = 26;
nbins = 30;

s = [58,82,88,149,190,595];

dmap = [[219,238,185]/255;...
    [151,213,187]/255;...
    [82,188,194]/255;...
    [75,150,189]/255;...
    [55,99,166]/255;...
    [41,54,144]/255];

cmap = lines(3);

burnin = 2e4;

MCMCmode = 'hierachichal';
%MCMCmode = 'single';

switch MCMCmode
    case 'hierachichal'
        % all genomic separation together
        vs = 1:6;
        load('Hierachichal_MCMC.mat')
        X = [X(:,1:5),X(:,5),X(:,6:end)];
    case 'single'
        % single genomic separation
        vs = 4;
        load(['Single',num2str(s(vs)),'_MCMC.mat']);
        X = [X(:,1:3),zeros(length(X),1),X(:,4),X(:,4),zeros(length(X),3),X(:,5:end)];
end

Xs = X(burnin:end,:);

m1X = mean(Xs,1);
s1X = std(Xs,1,1);
mdX = median(Xs,1);
p16X = prctile(Xs,16);
p84X = prctile(Xs,84);

%%% Plot Posterior distributions

H1=figure(1);
set(H1,'position',[50 700 2.2778*Wi 2*Le],'paperpositionmode','auto','color','w');
h11 = subplot(2,2,1,'parent',H1);
h12 = subplot(2,2,2,'parent',H1);
h13 = subplot(2,2,3,'parent',H1);
h14 = subplot(2,2,4,'parent',H1);
hold(h11,'on')
hold(h12,'on')
hold(h13,'on')
hold(h14,'on')
set(h11,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
set(h12,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
set(h13,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
set(h14,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
box(h11,'on')
box(h12,'on')
box(h13,'on')
box(h14,'on')

%b1
h=histogram(h11,Xs(:,1),nbins,'Normalization','probability');
h.FaceColor = cmap(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h11,'b_1 (1/min)')
ylabel(h11,'probability')
ylim(h11,[0,1.1*max(h.Values)])
title(h11,['b_1=',num2str(m1X(1),'%.3f'),'\pm',num2str(s1X(1),'%.3f')],'fontsize',ftsz)
%f2
h=histogram(h12,Xs(:,2),nbins,'Normalization','probability');
h.FaceColor = cmap(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h12,'f_2 (1/min)')
ylabel(h12,'probability')
ylim(h12,[0,1.1*max(h.Values)])
title(h12,['f_2=',num2str(m1X(2),'%.3f'),'\pm',num2str(s1X(2),'%.3f')],'fontsize',ftsz)
%b2
h=histogram(h13,Xs(:,3),nbins,'Normalization','probability');
h.FaceColor = cmap(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h13,'b_2 (1/min)')
ylabel(h13,'probability')
ylim(h13,[0,1.1*max(h.Values)])
title(h13,['b_2=',num2str(m1X(3),'%.3f'),'\pm',num2str(s1X(3),'%.3f')],'fontsize',ftsz)
%b3
h=histogram(h14,Xs(:,4),nbins,'Normalization','probability');
h.FaceColor = cmap(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h14,'b_3 (1/min)')
ylabel(h14,'probability')
ylim(h14,[0,1.1*max(h.Values)])
title(h14,['b_3=',num2str(m1X(4),'%.3f'),'\pm',num2str(s1X(4),'%.3f')],'fontsize',ftsz)

H2=figure(2);
set(H2,'position',[50 700 2.2778*Wi 2*Le],'paperpositionmode','auto','color','w');
h21 = subplot(2,2,1,'parent',H2);
h22 = subplot(2,2,2,'parent',H2);
h23 = subplot(2,2,3,'parent',H2);
h24 = subplot(2,2,4,'parent',H2);
hold(h21,'on')
hold(h22,'on')
hold(h23,'on')
hold(h24,'on')
set(h21,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
set(h22,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
set(h23,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
set(h24,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
box(h21,'on')
box(h22,'on')
box(h23,'on')
box(h24,'on')

%PD
h=histogram(h21,Xs(:,5),nbins,'Normalization','probability');
h.FaceColor = cmap(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h21,'\sigma_R(P) (nm)')
ylabel(h21,'probability')
ylim(h21,[0,1.1*max(h.Values)])
title(h21,['\sigma_R(P)=',num2str(m1X(5),'%.0f'),'\pm',num2str(s1X(5),'%.0f')],'fontsize',ftsz)
%red1
h=histogram(h22,Xs(:,7),nbins,'Normalization','probability');
h.FaceColor = cmap(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h22,'\sigma_I(O_{off})')
ylabel(h22,'probability')
ylim(h22,[0,1.1*max(h.Values)])
title(h22,['\sigma_I(O_{off})=',num2str(m1X(7),'%.3f'),'\pm',num2str(s1X(7),'%.3f')],'fontsize',ftsz)
%red2
h=histogram(h23,Xs(:,8),nbins,'Normalization','probability');
h.FaceColor = cmap(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h23,'\sigma_I(P_{off})')
ylabel(h23,'probability')
ylim(h23,[0,1.1*max(h.Values)])
title(h23,['\sigma_I(P_{off})=',num2str(m1X(8),'%.3f'),'\pm',num2str(s1X(8),'%.3f')],'fontsize',ftsz)
%red3
h=histogram(h24,Xs(:,9),nbins,'Normalization','probability');
h.FaceColor = cmap(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h24,'\sigma_I(P_{on})')
ylabel(h24,'probability')
ylim(h24,[0,1.1*max(h.Values)])
title(h24,['\sigma_I(P_{on})=',num2str(m1X(9),'%.3f'),'\pm',num2str(s1X(9),'%.3f')],'fontsize',ftsz)

for i=1:length(vs)
    Hi=figure(3+i);
    set(Hi,'position',[50 700 1*Wi 2*Le],'paperpositionmode','auto','color','w');
    hi1 = subplot(2,1,1,'parent',Hi);
    hi2 = subplot(2,1,2,'parent',Hi);
    hold(hi1,'on')
    hold(hi2,'on')
    set(hi1,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
    set(hi2,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
    box(hi1,'on')
    box(hi2,'on')
    
    %f1
    h=histogram(hi1,Xs(:,10+2*(i-1)),nbins,'Normalization','probability');
    h.FaceColor = dmap(i,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    ylim(hi1,[0,1.1*max(h.Values)])
    title(hi1,['f_1=',num2str(m1X(10+2*(i-1)),'%.3f'),'\pm',num2str(s1X(10+2*(i-1)),'%.3f')],'fontsize',ftsz)
    %OpenD
    h=histogram(hi2,Xs(:,11+2*(i-1)),nbins,'Normalization','probability');
    h.FaceColor = dmap(i,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    ylim(hi2,[0,1.1*max(h.Values)])
    title(hi2,['\sigma_R(O_{off})=',num2str(m1X(11+2*(i-1)),'%.0f'),'\pm',num2str(s1X(11+2*(i-1)),'%.0f')],'fontsize',ftsz)
    
    xlabel(hi1,'f_1 (1/min)')
    ylabel(hi1,'probability')
    
    xlabel(hi2,'\sigma_R(O_{off}) (nm)')
    ylabel(hi2,'probability')
end

H40=figure(40);
set(H40,'position',[50 700 2*Wi 1.5*Le],'paperpositionmode','auto','color','w');
h40 = subplot(1,1,1,'parent',H40);
hold(h40,'on')
set(h40,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
box(h40,'on')

h=histogram(h40,L(burnin:end,1),nbins,'Normalization','probability');
h.FaceColor = cmap(1,:);
h.EdgeColor = 'none';
h.FaceAlpha = 1;
xlabel(h40,'LogLikelihood')
ylabel(h40,'probability')
ylim(h40,[0,1.1*max(h.Values)])

%% Show MCMC trajectories
clc
close all

iteration = 1:length(L);
iteration = iteration(burnin:end);

H50=figure(50);
set(H50,'position',[50 700 Wi Le],'paperpositionmode','auto','color','w');
h50 = subplot(1,1,1,'parent',H50);
hold(h50,'on')
set(h50,'fontsize',ftsz,'linewidth',2)

plot(h50,L(:,1),'-','linewidth',1.5)
xlabel(h50,'iteration')
ylabel(h50,'LogLikelihood')

H1=figure(1);
set(H1,'position',[50 700 2*Wi 2*Le],'paperpositionmode','auto','color','w');
h11 = subplot(2,2,1,'parent',H1);
h12 = subplot(2,2,2,'parent',H1);
h13 = subplot(2,2,3,'parent',H1);
h14 = subplot(2,2,4,'parent',H1);
hold(h11,'on')
hold(h12,'on')
hold(h13,'on')
hold(h14,'on')
set(h11,'fontsize',ftsz,'linewidth',2)
set(h12,'fontsize',ftsz,'linewidth',2)
set(h13,'fontsize',ftsz,'linewidth',2)
set(h14,'fontsize',ftsz,'linewidth',2)

%b1
plot(h11,iteration,Xs(:,1),'-','linewidth',1.5)
xlabel(h11,'iteration')
ylabel(h11,'b_1 (1/min)')
%f2
plot(h12,iteration,Xs(:,2),'-','linewidth',1.5)
xlabel(h12,'iteration')
ylabel(h12,'f_2 (1/min)')
%b2
plot(h13,iteration,Xs(:,3),'-','linewidth',1.5)
xlabel(h13,'iteration')
ylabel(h13,'b_2 (1/min)')
%b3
plot(h14,iteration,Xs(:,4),'-','linewidth',1.5)
xlabel(h14,'iteration')
ylabel(h14,'b_3 (1/min)')

H2=figure(2);
set(H2,'position',[50 700 3*Wi 1*Le],'paperpositionmode','auto','color','w');
h21 = subplot(1,3,1,'parent',H2);
h22 = subplot(1,3,2,'parent',H2);
h23 = subplot(1,3,3,'parent',H2);
hold(h21,'on')
hold(h22,'on')
hold(h23,'on')
set(h21,'fontsize',ftsz,'linewidth',2)
set(h22,'fontsize',ftsz,'linewidth',2)
set(h23,'fontsize',ftsz,'linewidth',2)

%red1
plot(h21,iteration,Xs(:,7),'-','linewidth',1.5)
xlabel(h21,'iteration')
ylabel(h21,'\sigma_I(O_{off}) (a.u.)')
%red2
plot(h22,iteration,Xs(:,8),'-','linewidth',1.5)
xlabel(h22,'iteration')
ylabel(h22,'\sigma_I(P_{off}) (a.u.)')
%red3
plot(h23,iteration,Xs(:,9),'-','linewidth',1.5)
xlabel(h23,'iteration')
ylabel(h23,'\sigma_I(P_{on}) (a.u.)')

H3=figure(3);
set(H3,'position',[50 700 2*Wi 1*Le],'paperpositionmode','auto','color','w');
h31 = subplot(1,2,1,'parent',H3);
h32 = subplot(1,2,2,'parent',H3);
hold(h31,'on')
hold(h32,'on')
set(h31,'fontsize',ftsz,'linewidth',2)
set(h32,'fontsize',ftsz,'linewidth',2)

%PoffD
plot(h31,iteration,Xs(:,5),'-','linewidth',1.5)
xlabel(h31,'iteration')
ylabel(h31,'\sigma_R(P_{off}) (nm)')
%PonD
plot(h32,iteration,Xs(:,6),'-','linewidth',1.5)
xlabel(h32,'iteration')
ylabel(h32,'\sigma_R(P_{on}) (nm)')

for i=1:length(vs)
    Hi=figure(3+i);
    set(Hi,'position',[50 700 2*Wi 1*Le],'paperpositionmode','auto','color','w');
    sgtitle(Hi,['s=',num2str(s(vs(i))),' kb'],'fontsize',ftsz,'fontweight','bold')
    hi1 = subplot(1,2,1,'parent',Hi);
    hi2 = subplot(1,2,2,'parent',Hi);
    hold(hi1,'on')
    hold(hi2,'on')
    set(hi1,'fontsize',ftsz,'linewidth',2)
    set(hi2,'fontsize',ftsz,'linewidth',2)
    
    %f1
    plot(hi1,iteration,Xs(:,10+2*(i-1)),'-','linewidth',1.5)
    xlabel(hi1,'iteration')
    ylabel(hi1,'f_1 (1/min)')
    %OpenD
    plot(hi2,iteration,Xs(:,11+2*(i-1)),'-','linewidth',1.5)
    xlabel(hi2,'iteration')
    ylabel(hi2,'\sigma_R(O_{off}) (nm)')
end
