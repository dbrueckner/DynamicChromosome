%CheckObservation
%   Plot posterior distributions from MCMC samples
%
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

clear
clc
close all

pobsR = @(r,na) (3/(2*pi*na^2))^(3/2)*4*pi*r.^2.*exp(-3*r.^2/(2*na^2));
pobsI = @(x,b) x./b^2 .*exp(-x.^2./(2.*b^2));

Wi = 360;
Le = 330;

ftsz = 26;
lnwi = 5;

timeRes = 0.46;

load('fullDataSet.mat')

burnin = 2e4;
Xs = results(burnin:end,:);
m1X = mean(Xs,1);

vs = 1:6;
s = [57.8,81.7,87.5,148.7,189.6,595.1];

xr1 = linspace(0,3e3,100);
xr2 = linspace(0,1.2e3,100);
xr3 = linspace(0,1.2e3,100);

i1edges = linspace(0,0.6,21);
i2edges = linspace(0,0.6,21);
i3edges = linspace(0,5,31);
xi1 = linspace(0,1,100);
xi2 = linspace(0,1,100);
xi3 = linspace(0,5,100);

smap = [[0.5,0.5,0.5];...
    [0,0,1];...
    [0,1,1];...
    [1,0,0]];

dmap = [[219,238,185]/255;...
    [151,213,187]/255;...
    [82,188,194]/255;...
    [75,150,189]/255;...
    [55,99,166]/255;...
    [41,54,144]/255];

for k=vs
    smatrix = zeros(length(state),99);
    dmatrix = zeros(length(state),99);
    imatrix = zeros(length(state),99);
    for i=1:length(state)
        st = state{k,i};
        dist = bgDist{k,i};
        inte = red{k,i};
        smatrix(i,1:length(st)) = 1+st;
        dmatrix(i,1:length(st)) = dist;
        imatrix(i,1:length(st)) = inte;
    end
    
    I1 = smatrix==1;
    I2 = smatrix==2;
    I3 = smatrix==3;
    
    H1=figure(1+3*(k-1));
    set(H1,'position',[50 700 Wi 3*Le],'paperpositionmode','auto','color','w');
    %sgtitle(H1,['s=',num2str(s(k)),' kb'],'fontsize',ftsz,'fontweight','bold')
    h11 = subplot(3,1,1,'parent',H1);
    h12 = subplot(3,1,2,'parent',H1);
    h13 = subplot(3,1,3,'parent',H1);
    hold(h11,'on')
    hold(h12,'on')
    hold(h13,'on')
    set(h11,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
    set(h12,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
    set(h13,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
    box(h11,'on')
    box(h12,'on')
    box(h13,'on')
    xlim(h11,[0,2.5e3])
    xlim(h12,[0,1.2e3])
    xlim(h13,[0,1.2e3])
    
    h=histogram(h11,dmatrix(I1),'Normalization','pdf','BinMethod','auto');
    h.FaceColor = dmap(k,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    %plot(h11,h.BinEdges(1:end-1)+h.BinWidth/2,h.Values,'-','linewidth',3,'color',dmap(k,:))
    plot(h11,xr1,pobsR(xr1,m1X(11+2*(k-1))),':','linewidth',lnwi,'color',smap(2,:))
    xlabel(h11,'R (nm)')
    ylabel(h11,'pdf')
    ylim(h11,[0,1.1*max(max(h.Values),max(pobsR(xr1,m1X(11+2*(k-1)))))])
    title(h11,'Distance in O_{off}','fontsize',ftsz)
    
    h=histogram(h12,dmatrix(I2),'Normalization','pdf','BinMethod','auto');
    h.FaceColor = dmap(k,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    %plot(h12,h.BinEdges(1:end-1)+h.BinWidth/2,h.Values,'-','linewidth',3,'color',dmap(k,:))
    plot(h12,xr2,pobsR(xr2,m1X(5)),':','linewidth',lnwi,'color',smap(3,:))
    xlabel(h12,'R (nm)')
    ylabel(h12,'pdf')
    ymax = max(max(h.Values),max(pobsR(xr2,m1X(5))));
    ymax(isnan(ymax)) = 1;
    ylim(h12,[0,1.1*ymax])
    title(h12,'Distance in P_{off}','fontsize',ftsz)
    
    h=histogram(h13,dmatrix(I3),'Normalization','pdf','BinMethod','auto');
    h.FaceColor = dmap(k,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    %plot(h13,h.BinEdges(1:end-1)+h.BinWidth/2,h.Values,'-','linewidth',3,'color',dmap(k,:))
    plot(h13,xr3,pobsR(xr3,m1X(6)),':','linewidth',lnwi,'color',smap(4,:))
    xlabel(h13,'R (nm)')
    ylabel(h13,'pdf')
    ymax = max(max(h.Values),max(pobsR(xr3,m1X(6))));
    ymax(isnan(ymax)) = 1;
    ylim(h13,[0,1.1*ymax])
    title(h13,'Distance in P_{on}','fontsize',ftsz)
    
    H2=figure(2+3*(k-1));
    set(H2,'position',[50 700 Wi 3*Le],'paperpositionmode','auto','color','w');
    %sgtitle(H2,['s=',num2str(s(k)),' kb'],'fontsize',ftsz,'fontweight','bold')
    h21 = subplot(3,1,1,'parent',H2);
    h22 = subplot(3,1,2,'parent',H2);
    h23 = subplot(3,1,3,'parent',H2);
    hold(h21,'on')
    hold(h22,'on')
    hold(h23,'on')
    set(h21,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
    set(h22,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
    set(h23,'fontsize',ftsz,'linewidth',2,'tickdir','out','ytick',[])
    box(h21,'on')
    box(h22,'on')
    box(h23,'on')
    xlim(h21,[0,0.6])
    xlim(h22,[0,0.6])
    xlim(h23,[0,4])
    
    h=histogram(h21,imatrix(I1),i1edges,'Normalization','pdf');
    h.FaceColor = dmap(k,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    plot(h21,xi1,pobsI(xi1,m1X(7)),':','linewidth',lnwi,'color',smap(2,:))
    xlabel(h21,'I (a.u.)')
    ylabel(h21,'pdf')
    ylim(h21,[0,1.1*max(max(h.Values),max(pobsI(xi1,m1X(7))))])
    title(h21,'Intensity in O_{off}','fontsize',ftsz)
    
    h=histogram(h22,imatrix(I2),i2edges,'Normalization','pdf');
    h.FaceColor = dmap(k,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    plot(h22,xi2,pobsI(xi2,m1X(8)),':','linewidth',lnwi,'color',smap(3,:))
    xlabel(h22,'I (a.u.)')
    ylabel(h22,'pdf')
    ymax = max(max(h.Values),max(pobsI(xi2,m1X(8))));
    ymax(isnan(ymax)) = 1;
    ylim(h22,[0,1.1*ymax])
    title(h22,'Intensity in P_{off}','fontsize',ftsz)
    
    h=histogram(h23,imatrix(I3),i3edges,'Normalization','pdf');
    h.FaceColor = dmap(k,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    plot(h23,xi3,pobsI(xi3,m1X(9)),':','linewidth',lnwi,'color',smap(4,:))
    xlabel(h23,'I (a.u.)')
    ylabel(h23,'pdf')
    ymax = max(max(h.Values),max(pobsI(xi3,m1X(9))));
    ymax(isnan(ymax)) = 1;
    ylim(h23,[0,1.1*ymax])
    title(h23,'Intensity in P_{on}','fontsize',ftsz)
    
    H3=figure(3+3*(k-1));
    set(H3,'position',[50 700 1.4*Wi 2*Le],'paperpositionmode','auto','color','w');
    h3 = subplot(1,1,1,'parent',H3);
    hold(h3,'on')
    
    set(h3,'fontsize',22,'linewidth',2,'tickdir','out','ydir','reverse','layer','top')
    xlabel(h3,'time (min)')
    ylabel(h3,'cell index')
    
    j=find(smatrix(:,1)==0,1,'first');
    if isempty(j)
        j=length(state);
    else
        j=j-1;
    end
    
    sum(I3(:))/(sum(I1(:))+sum(I2(:))+sum(I3(:)))
    
    imagesc(h3,(0:98)*timeRes,1:j,smatrix(1:j,:))
    colormap(h3,smap)
    xlim(h3,[0,35])
    ylim(h3,[1,j])
    box(h3,'on')
end