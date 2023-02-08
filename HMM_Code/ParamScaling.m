%ParamScaling
%   Investigate the HMM parameter scaling
%
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

clear
close all
clc

Wi = 500;
Le = 375;

ftsz = 26;
mksz = 12;
lnwi = 3;

burnin = 2e4;
load('Hierachichal_MCMC.mat','X')
X = [X(:,1:5),X(:,5),X(:,6:end)];
X_global = X(burnin:end,:);
mdX_global = median(X_global);
plX_global = prctile(X_global,16);
puX_global = prctile(X_global,84);
vlogX_global = var(log(X_global));

mylist = {'58','82','88','149','190','595'};
mydist = [57.8,81.7,87.5,148.7,189.6,595.1];
cmap = [[219,238,185]/255;...
    [151,213,187]/255;...
    [82,188,194]/255;...
    [75,150,189]/255;...
    [55,99,166]/255;...
    [41,54,144]/255];

burnin = 1e4;

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
set(h11,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
set(h12,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
set(h13,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
set(h14,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h11,'s (kb)')
xlabel(h12,'s (kb)')
xlabel(h13,'s (kb)')
xlabel(h14,'s (kb)')
ylabel(h11,'b_1 (1/min)')
ylabel(h12,'f_2 (1/min)')
ylabel(h13,'b_2 (1/min)')
ylabel(h14,'\sigma_R(P) (nm)')
xlim(h11,[3e1,1e3])
xlim(h12,[3e1,1e3])
xlim(h13,[3e1,1e3])
xlim(h14,[3e1,1e3])
box(h11,'on')
box(h12,'on')
box(h13,'on')
box(h14,'on')

plot(h11,[3e1,1e3],mdX_global(1)*[1,1],'-k','linewidth',lnwi)
plot(h11,[3e1,1e3],plX_global(1)*[1,1],':k','linewidth',lnwi)
plot(h11,[3e1,1e3],puX_global(1)*[1,1],':k','linewidth',lnwi)

plot(h12,[3e1,1e3],mdX_global(2)*[1,1],'-k','linewidth',lnwi)
plot(h12,[3e1,1e3],plX_global(2)*[1,1],':k','linewidth',lnwi)
plot(h12,[3e1,1e3],puX_global(2)*[1,1],':k','linewidth',lnwi)

plot(h13,[3e1,1e3],mdX_global(3)*[1,1],'-k','linewidth',lnwi)
plot(h13,[3e1,1e3],plX_global(3)*[1,1],':k','linewidth',lnwi)
plot(h13,[3e1,1e3],puX_global(3)*[1,1],':k','linewidth',lnwi)

plot(h14,[3e1,1e3],mdX_global(5)*[1,1],'-k','linewidth',lnwi)
plot(h14,[3e1,1e3],plX_global(5)*[1,1],':k','linewidth',lnwi)
plot(h14,[3e1,1e3],puX_global(5)*[1,1],':k','linewidth',lnwi)

H2=figure(2);
set(H2,'position',[50 700 2*Wi 0.9280*Le],'paperpositionmode','auto','color','w');
h21 = subplot(1,2,1,'parent',H2);
h22 = subplot(1,2,2,'parent',H2);
hold(h21,'on')
hold(h22,'on')
set(h21,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
set(h22,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h21,'s (kb)')
xlabel(h22,'s (kb)')
xlim(h21,[3e1,1e3])
xlim(h22,[3e1,1e3])
ylabel(h21,'f_1 (1/min)')
ylabel(h22,'\sigma_R(O_{off}) (nm)')
box(h21,'on')
box(h22,'on')

H3=figure(3);
set(H3,'position',[50 700 2*Wi 2*Le],'paperpositionmode','auto','color','w');
h31 = subplot(2,2,1,'parent',H3);
h32 = subplot(2,2,2,'parent',H3);
h33 = subplot(2,2,4,'parent',H3);
h34 = subplot(2,2,3,'parent',H3);
hold(h31,'on')
hold(h32,'on')
hold(h33,'on')
hold(h34,'on')
set(h31,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
set(h32,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
set(h33,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
set(h34,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h31,'s (kb)')
xlabel(h32,'s (kb)')
xlabel(h33,'s (kb)')
xlabel(h34,'s (kb)')
ylabel(h31,'P(O_{off})')
ylabel(h32,'P(P_{off})')
ylabel(h33,'P(P_{on})')
ylabel(h34,'P(P)')
xlim(h31,[3e1,1e3])
xlim(h32,[3e1,1e3])
xlim(h33,[3e1,1e3])
xlim(h34,[3e1,1e3])
ylim(h31,[0,1])
ylim(h32,[0,1])
ylim(h33,[0,1])
box(h31,'on')
box(h32,'on')
box(h33,'on')
box(h34,'on')

H4=figure(4);
set(H4,'position',[50 700 2*Wi 0.9280*Le],'paperpositionmode','auto','color','w');
h41 = subplot(1,2,1,'parent',H4);
h42 = subplot(1,2,2,'parent',H4);
hold(h41,'on')
hold(h42,'on')
set(h41,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
set(h42,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h41,'s (kb)')
xlabel(h42,'s (kb)')
ylabel(h41,'P(P_{off}|P)')
ylabel(h42,'P(P_{on}|P)')
xlim(h41,[3e1,1e3])
xlim(h42,[3e1,1e3])
ylim(h41,[0,1])
ylim(h42,[0,1])
box(h41,'on')
box(h42,'on')

Nd = length(mylist);

X = log(mydist);

Y11 = nan(Nd,1);
V11 = nan(Nd,1);

Y12 = nan(Nd,1);
V12 = nan(Nd,1);

Y13 = nan(Nd,1);
V13 = nan(Nd,1);

Y14 = nan(Nd,1);
V14 = nan(Nd,1);

Y21 = nan(Nd,1);
V21 = nan(Nd,1);
Y21_global = nan(Nd,1);
V21_global = nan(Nd,1);

Y22 = nan(Nd,1);
V22 = nan(Nd,1);
Y22_global = nan(Nd,1);
V22_global = nan(Nd,1);

Y31 = nan(Nd,1);
V31 = nan(Nd,1);
EB31 = nan(Nd,2);
Y31_global = nan(Nd,1);
V31_global = nan(Nd,1);
EB31_global = nan(Nd,2);

Y32 = nan(Nd,1);
V32 = nan(Nd,1);
EB32 = nan(Nd,2);
Y32_global = nan(Nd,1);
V32_global = nan(Nd,1);
EB32_global = nan(Nd,2);

Y33 = nan(Nd,1);
V33 = nan(Nd,1);
EB33 = nan(Nd,2);
Y33_global = nan(Nd,1);
V33_global = nan(Nd,1);
EB33_global = nan(Nd,2);

Y34 = nan(Nd,1);
V34 = nan(Nd,1);
Y34_global = nan(Nd,1);
V34_global = nan(Nd,1);

for i=1:Nd
    load(['Single',mylist{i},'_MCMC'],'X')
    X_local = X(burnin:end,:);
    
    mdX = median(X_local);
    plX = prctile(X_local,16);
    puX = prctile(X_local,84);
    logvX = var(log(X_local));
    
    P_global = GetOccupancy([X_global(:,10+2*(i-1)),X_global(:,1:4)]);
    P_global = [P_global,P_global(:,2)+P_global(:,3)];
    Pl_global = prctile(P_global,16);
    Pu_global = prctile(P_global,84);
    Pm_global = median(P_global);
    logvP_global = var(log(P_global));
    
    P_local = GetOccupancy([X_local(:,5),X_local(:,1:3),zeros(size(X_local,1),1)]);
    P_local = [P_local,P_local(:,2)+P_local(:,3)];
    Pl_local = prctile(P_local,16);
    Pu_local = prctile(P_local,84);
    Pm_local = median(P_local);
    logvP_local = var(log(P_local));
    
    %b1
    errorbar(h11,mydist(i),mdX(1),mdX(1)-plX(1),puX(1)-mdX(1),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y11(i) = log(mdX(1));
    V11(i) = logvX(1);
    %f2
    errorbar(h12,mydist(i),mdX(2),mdX(2)-plX(2),puX(2)-mdX(2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y12(i) = log(mdX(2));
    V12(i) = logvX(2);
    %b2
    errorbar(h13,mydist(i),mdX(3),mdX(3)-plX(3),puX(3)-mdX(3),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y13(i) = log(mdX(3));
    V13(i) = logvX(3);
    %PD
    errorbar(h14,mydist(i),mdX(4),mdX(4)-plX(4),puX(4)-mdX(4),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y14(i) = log(mdX(4));
    V14(i) = logvX(4);
    
    %f1
    errorbar(h21,mydist(i),mdX_global(10+2*(i-1)),mdX_global(10+2*(i-1))-plX_global(10+2*(i-1)),puX_global(10+2*(i-1))-mdX_global(10+2*(i-1)),...
        'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h21,mydist(i),mdX(5),mdX(5)-plX(5),puX(5)-mdX(5),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y21(i) = log(mdX(5));
    V21(i) = logvX(5);
    Y21_global(i) = log(mdX_global(10+2*(i-1)));
    V21_global(i) = vlogX_global(10+2*(i-1));
    
    %OpenD
    errorbar(h22,mydist(i),mdX_global(11+2*(i-1)),mdX_global(11+2*(i-1))-plX_global(11+2*(i-1)),puX_global(11+2*(i-1))-mdX_global(11+2*(i-1)),...
        'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h22,mydist(i),mdX(6),mdX(6)-plX(6),puX(6)-mdX(6),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y22(i) = log(mdX(6));
    V22(i) = logvX(6);
    Y22_global(i) = log(mdX_global(11+2*(i-1)));
    V22_global(i) = vlogX_global(11+2*(i-1));
    
    %P(open)
    errorbar(h31,mydist(i),Pm_global(1),Pm_global(1)-Pl_global(1),Pu_global(1)-Pm_global(1),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h31,mydist(i),Pm_local(1),Pm_local(1)-Pl_local(1),Pu_local(1)-Pm_local(1),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y31(i) = log(Pm_local(1));
    V31(i) = logvP_local(1);
    EB31(i,:) = [Pm_local(1)-Pl_local(1),Pu_local(1)-Pm_local(1)];
    Y31_global(i) = log(Pm_global(1));
    V31_global(i) = logvP_global(1);
    EB31_global(i,:) = [Pm_global(1)-Pl_global(1),Pu_global(1)-Pm_global(1)];
    
    %P(off)
    errorbar(h32,mydist(i),Pm_global(2),Pm_global(2)-Pl_global(2),Pu_global(2)-Pm_global(2),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h32,mydist(i),Pm_local(2),Pm_local(2)-Pl_local(2),Pu_local(2)-Pm_local(2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y32(i) = log(Pm_local(2));
    V32(i) = logvP_local(2);
    EB32(i,:) = [Pm_local(2)-Pl_local(2),Pu_local(2)-Pm_local(2)];
    Y32_global(i) = log(Pm_global(2));
    V32_global(i) = logvP_global(2);
    EB32_global(i,:) = [Pm_global(2)-Pl_global(2),Pu_global(2)-Pm_global(2)];
    
    %P(on)
    errorbar(h33,mydist(i),Pm_global(3),Pm_global(3)-Pl_global(3),Pu_global(3)-Pm_global(3),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h33,mydist(i),Pm_local(3),Pm_local(3)-Pl_local(3),Pu_local(3)-Pm_local(3),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y33(i) = log(Pm_local(3));
    V33(i) = logvP_local(3);
    EB33(i,:) = [Pm_local(3)-Pl_local(3),Pu_local(3)-Pm_local(3)];
    Y33_global(i) = log(Pm_global(3));
    V33_global(i) = logvP_global(3);
    EB33_global(i,:) = [Pm_global(3)-Pl_global(3),Pu_global(3)-Pm_global(3)];
    
    %P(P)
    errorbar(h34,mydist(i),Pm_global(4),Pm_global(4)-Pl_global(4),Pu_global(4)-Pm_global(4),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h34,mydist(i),Pm_local(4),Pm_local(4)-Pl_local(4),Pu_local(4)-Pm_local(4),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y34(i) = log(Pm_local(4));
    V34(i) = logvP_local(4);
    Y34_global(i) = log(Pm_global(4));
    V34_global(i) = logvP_global(4);
    
    Pp_global = P_global(:,2:3)./repmat((P_global(:,2)+P_global(:,3)),1,2);
    Ppl_global = prctile(Pp_global,16);
    Ppu_global = prctile(Pp_global,84);
    Ppm_global = median(Pp_global);
    
    Pp_local = P_local(:,2:3)./repmat((P_local(:,2)+P_local(:,3)),1,2);
    Ppl_local = prctile(Pp_local,16);
    Ppu_local = prctile(Pp_local,84);
    Ppm_local = median(Pp_local);
    
    %P(off|P)
    errorbar(h41,mydist(i),Ppm_local(1),Ppm_local(1)-Ppl_local(1),Ppu_local(1)-Ppm_local(1),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    
    %P(on|P)
    errorbar(h42,mydist(i),Ppm_local(2),Ppm_local(2)-Ppl_local(2),Ppu_local(2)-Ppm_local(2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
end

p=plot(h41,[3e1,1e3],Ppm_global(1)*[1,1],'-k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h41,[3e1,1e3],Ppl_global(1)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h41,[3e1,1e3],Ppu_global(1)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')

p=plot(h42,[3e1,1e3],Ppm_global(2)*[1,1],'-k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h42,[3e1,1e3],Ppl_global(2)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h42,[3e1,1e3],Ppu_global(2)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')

y=ylim(h11);
ylim(h11,[0,y(2)])
y=ylim(h12);
ylim(h12,[0,y(2)])
y=ylim(h13);
ylim(h13,[0,y(2)])

x = logspace(log10(3e1),log10(1e3),100);

fopt1 = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[-Inf,-Inf],...
    'Upper',[+Inf,+Inf],...
    'Startpoint',[0,0]);

ftype1 = fittype('a*x+b','options',fopt1);

ft = fit(log(mydist(:)),Y11(:),ftype1,'Weight',1./V11(:));
CI11 = confint(ft,0.68);
P11 = [ft.a,ft.b];
E11 = mean(abs(CI11-repmat(P11,2,1)));
p=plot(h11,x,exp(P11(2))*x.^P11(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y12(:),ftype1,'Weight',1./V12(:));
CI12 = confint(ft,0.68);
P12 = [ft.a,ft.b];
E12 = mean(abs(CI12-repmat(P12,2,1)));
p=plot(h12,x,exp(P12(2))*x.^P12(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y13(:),ftype1,'Weight',1./V13(:));
CI13 = confint(ft,0.68);
P13 = [ft.a,ft.b];
E13 = mean(abs(CI13-repmat(P13,2,1)));
p=plot(h13,x,exp(P13(2))*x.^P13(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y14(:),ftype1,'Weight',1./V14(:));
CI14 = confint(ft,0.68);
P14 = [ft.a,ft.b];
E14 = mean(abs(CI14-repmat(P14,2,1)));
p=plot(h14,x,exp(P14(2))*x.^P14(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

annotation(H1,'textbox',...
    [0.33 0.83 0.067 0.085],...
    'String',['~s^{',num2str(P11(1),'%.2f'),'\pm',num2str(E11(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H1,'textbox',...
    [0.77 0.83 0.067 0.085],...
    'String',['~s^{',num2str(P12(1),'%.2f'),'\pm',num2str(E12(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H1,'textbox',...
    [0.33 0.36 0.067 0.085],...
    'String',['~s^{',num2str(P13(1),'%.2f'),'\pm',num2str(E13(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H1,'textbox',...
    [0.77 0.36 0.067 0.085],...
    'String',['~s^{',num2str(P14(1),'%.2f'),'\pm',num2str(E14(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

ft = fit(log(mydist(:)),Y21(:),ftype1,'Weight',1./V21(:));
CI21 = confint(ft,0.68);
P21 = [ft.a,ft.b];
E21 = mean(abs(CI21-repmat(P21,2,1)));
p=plot(h21,x,exp(P21(2))*x.^P21(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

xx = log(mydist(1:(end-1)));
yy = Y22(1:(end-1));
ww = 1./V22(1:(end-1));
ft = fit(xx(:),yy(:),ftype1,'Weight',ww(:));
%ft = fit(log(mydist(:)),Y22(:),ftype1,'Weight',1./V22(:));
CI22 = confint(ft,0.68);
P22 = [ft.a,ft.b];
E22 = mean(abs(CI22-repmat(P22,2,1)));
p=plot(h22,x,exp(P22(2))*x.^P22(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y21_global(:),ftype1,'Weight',1./V21_global(:));
CI21_global = confint(ft,0.68);
P21_global = [ft.a,ft.b];
E21_global = mean(abs(CI21_global-repmat(P21_global,2,1)));
p=plot(h21,x,exp(P21_global(2))*x.^P21_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
uistack(p,'bottom')

xx = log(mydist(1:(end-1)));
yy = Y22_global(1:(end-1));
ww = 1./V22_global(1:(end-1));
ft = fit(xx(:),yy(:),ftype1,'Weight',ww(:));
CI22_global = confint(ft,0.68);
P22_global = [ft.a,ft.b];
E22_global = mean(abs(CI22_global-repmat(P22_global,2,1)));
p=plot(h22,x,exp(P22_global(2))*x.^P22_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
uistack(p,'bottom')

annotation(H2,'textbox',...
    [0.33 0.82 0.067 0.085],...
    'String',['~s^{',num2str(P21_global(1),'%.2f'),'\pm',num2str(E21_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H2,'textbox',...
    [0.59 0.82 0.067 0.085],...
    'String',['~s^{',num2str(P22_global(1),'%.2f'),'\pm',num2str(E22_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

ft = fit(log(mydist(:)),Y32(:),ftype1,'Weight',1./V32(:));
CI32 = confint(ft,0.68);
P32 = [ft.a,ft.b];
E32 = mean(abs(CI32-repmat(P32,2,1)));
p=plot(h32,x,exp(P32(2))*x.^P32(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y33(:),ftype1,'Weight',1./V33(:));
CI33 = confint(ft,0.68);
P33 = [ft.a,ft.b];
E33 = mean(abs(CI33-repmat(P33,2,1)));
p=plot(h33,x,exp(P33(2))*x.^P33(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y34(:),ftype1,'Weight',1./V34(:));
CI34 = confint(ft,0.68);
P34 = [ft.a,ft.b];
E34 = mean(abs(CI34-repmat(P34,2,1)));
p=plot(h34,x,exp(P34(2))*x.^P34(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y32_global(:),ftype1,'Weight',1./V32_global(:));
CI32_global = confint(ft,0.68);
P32_global = [ft.a,ft.b];
E32_global = mean(abs(CI32_global-repmat(P32_global,2,1)));
p=plot(h32,x,exp(P32_global(2))*x.^P32_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y33_global(:),ftype1,'Weight',1./V33_global(:));
CI33_global = confint(ft,0.68);
P33_global = [ft.a,ft.b];
E33_global = mean(abs(CI33_global-repmat(P33_global,2,1)));
p=plot(h33,x,exp(P33_global(2))*x.^P33_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y34_global(:),ftype1,'Weight',1./V34_global(:));
CI34_global = confint(ft,0.68);
P34_global = [ft.a,ft.b];
E34_global = mean(abs(CI34_global-repmat(P34_global,2,1)));
p=plot(h34,x,exp(P34_global(2))*x.^P34_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
uistack(p,'bottom')

annotation(H3,'textbox',...
    [0.77 0.83 0.067 0.085],...
    'String',['~s^{',num2str(P32_global(1),'%.2f'),'\pm',num2str(E32_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H3,'textbox',...
    [0.77 0.36 0.067 0.085],...
    'String',['~s^{',num2str(P33_global(1),'%.2f'),'\pm',num2str(E33_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H3,'textbox',...
    [0.33 0.36 0.067 0.085],...
    'String',['~s^{',num2str(P34_global(1),'%.2f'),'\pm',num2str(E34_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% Lifetime

H5=figure(5);
set(H5,'position',[50 700 2*Wi 2*Le],'paperpositionmode','auto','color','w');
h51 = subplot(2,2,1,'parent',H5);
h52 = subplot(2,2,2,'parent',H5);
h53 = subplot(2,2,4,'parent',H5);
h54 = subplot(2,2,3,'parent',H5);
hold(h51,'on')
hold(h52,'on')
hold(h53,'on')
hold(h54,'on')
set(h51,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
set(h52,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
set(h53,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
set(h54,'fontsize',ftsz,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h51,'s (kb)')
xlabel(h52,'s (kb)')
xlabel(h53,'s (kb)')
xlabel(h54,'s (kb)')
ylabel(h51,'T(O_{off}) (min)')
ylabel(h52,'T(P_{off}) (min)')
ylabel(h53,'T(P_{on}) (min)')
ylabel(h54,'T(P) (min)')
xlim(h51,[3e1,1e3])
xlim(h52,[3e1,1e3])
xlim(h53,[3e1,1e3])
xlim(h54,[3e1,1e3])
box(h51,'on')
box(h52,'on')
box(h53,'on')
box(h54,'on')

H20=figure(20);
set(H20,'position',[50 700 Wi Le],'paperpositionmode','auto','color','w');
h20 = subplot(1,1,1,'parent',H20);
hold(h20,'on')
set(h20,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h20,'s (kb)')
ylabel(h20,'F(P_{on}) (1/h)')
xlim(h20,[3e1,1e3])
box(h20,'on')

Y51 = nan(Nd,1);
V51 = nan(Nd,1);
EB51 = nan(Nd,2);
Y51_global = nan(Nd,1);
V51_global = nan(Nd,1);
EB51_global = nan(Nd,2);

Y52 = nan(Nd,1);
V52 = nan(Nd,1);
EB52 = nan(Nd,2);

Y53 = nan(Nd,1);
V53 = nan(Nd,1);
EB53 = nan(Nd,2);

Y54 = nan(Nd,1);
V54 = nan(Nd,1);
EB54 = nan(Nd,2);

Y20_global = nan(Nd,1);
V20_global = nan(Nd,1);
Y20 = nan(Nd,1);
V20 = nan(Nd,1);

for i=1:Nd
    load(['Single',mylist{i},'_MCMC'],'X')
    X_local = X(burnin:end,:);
    
    T_global = GetLifetime([X_global(:,10+2*(i-1)),X_global(:,1:4)]);
    Tl_global = prctile(T_global,16);
    Tu_global = prctile(T_global,84);
    Tm_global = median(T_global);
    logvT_global = var(log(T_global));
    
    T_local = GetLifetime([X_local(:,5),X_local(:,1:3),zeros(size(X_local,1),1)]);
    Tl_local = prctile(T_local,16);
    Tu_local = prctile(T_local,84);
    Tm_local = median(T_local);
    logvT_local = var(log(T_local));
    
    errorbar(h51,mydist(i),Tm_global(1),Tm_global(1)-Tl_global(1),Tu_global(1)-Tm_global(1),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h51,mydist(i),Tm_local(1),Tm_local(1)-Tl_local(1),Tu_local(1)-Tm_local(1),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y51(i) = log(Tm_local(1));
    V51(i) = logvT_local(1);
    EB51(i,:) = [Tm_local(1)-Tl_local(1),Tu_local(1)-Tm_local(1)];
    Y51_global(i) = log(Tm_global(1));
    V51_global(i) = logvT_global(1);
    EB51_global(i,:) = [Tm_global(1)-Tl_global(1),Tu_global(1)-Tm_global(1)];
    
    errorbar(h52,mydist(i),Tm_local(2),Tm_local(2)-Tl_local(2),Tu_local(2)-Tm_local(2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y52(i) = log(Tm_local(2));
    V52(i) = logvT_local(2);
    EB52(i,:) = [Tm_local(2)-Tl_local(2),Tu_local(2)-Tm_local(2)];
    
    errorbar(h53,mydist(i),Tm_local(3),Tm_local(3)-Tl_local(3),Tu_local(3)-Tm_local(3),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y53(i) = log(Tm_local(3));
    V53(i) = logvT_local(3);
    EB53(i,:) = [Tm_local(3)-Tl_local(3),Tu_local(3)-Tm_local(3)];
    
    errorbar(h54,mydist(i),Tm_local(4),Tm_local(4)-Tl_local(4),Tu_local(4)-Tm_local(4),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y54(i) = log(Tm_local(4));
    V54(i) = logvT_local(4);
    EB54(i,:) = [Tm_local(4)-Tl_local(4),Tu_local(4)-Tm_local(4)];
    
    Fon_global = 60./(T_global(:,3)+T_global(:,5));
    Fonl_global = prctile(Fon_global,16);
    Fonu_global = prctile(Fon_global,84);
    Fonm_global = median(Fon_global);
    logvFon_global = var(log(Fon_global));
    
    Fon_local = 60./(T_local(:,3)+T_local(:,5));
    Fonl_local = prctile(Fon_local,16);
    Fonu_local = prctile(Fon_local,84);
    Fonm_local = median(Fon_local);
    logvFon_local = var(log(Fon_local));
    
    errorbar(h20,mydist(i),Fonm_global(1),Fonm_global(1)-Fonl_global(1),Fonu_global(1)-Fonm_global(1),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h20,mydist(i),Fonm_local(1),Fonm_local(1)-Fonl_local(1),Fonu_local(1)-Fonm_local(1),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y20(i) = log(Fonm_local(1));
    V20(i) = logvFon_local(1);
    Y20_global(i) = log(Fonm_global(1));
    V20_global(i) = logvFon_global(1);
    
end

p=plot(h52,[3e1,1e3],Tm_global(2)*[1,1],'-k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h52,[3e1,1e3],Tl_global(2)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h52,[3e1,1e3],Tu_global(2)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')
y=ylim(h52);
ylim(h52,[0,y(2)])

p=plot(h53,[3e1,1e3],Tm_global(3)*[1,1],'-k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h53,[3e1,1e3],Tl_global(3)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h53,[3e1,1e3],Tu_global(3)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')

p=plot(h54,[3e1,1e3],Tm_global(4)*[1,1],'-k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h54,[3e1,1e3],Tl_global(4)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h54,[3e1,1e3],Tu_global(4)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')

x = logspace(log10(3e1),log10(1e3),100);

ft = fit(log(mydist(:)),Y51(:),ftype1,'Weight',1./V51(:));
CI51 = confint(ft,0.68);
P51 = [ft.a,ft.b];
E51 = mean(abs(CI51-repmat(P51,2,1)));
p=plot(h51,x,exp(P51(2))*x.^P51(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y52(:),ftype1,'Weight',1./V52(:));
CI52 = confint(ft,0.68);
P52 = [ft.a,ft.b];
E52 = mean(abs(CI52-repmat(P52,2,1)));
p=plot(h52,x,exp(P52(2))*x.^P52(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y53(:),ftype1,'Weight',1./V53(:));
CI53 = confint(ft,0.68);
P53 = [ft.a,ft.b];
E53 = mean(abs(CI53-repmat(P53,2,1)));
p=plot(h53,x,exp(P53(2))*x.^P53(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y54(:),ftype1,'Weight',1./V54(:));
CI54 = confint(ft,0.68);
P54 = [ft.a,ft.b];
E54 = mean(abs(CI54-repmat(P54,2,1)));
p=plot(h54,x,exp(P54(2))*x.^P54(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y51_global(:),ftype1,'Weight',1./V51_global(:));
CI51_global = confint(ft,0.68);
P51_global = [ft.a,ft.b];
E51_global = mean(abs(CI51_global-repmat(P51_global,2,1)));
p=plot(h51,x,exp(P51_global(2))*x.^P51_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
uistack(p,'bottom')

annotation(H5,'textbox',...
    [0.15 0.835 0.067 0.085],...
    'String',['~s^{',num2str(P51_global(1),'%.2f'),'\pm',num2str(E51_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H5,'textbox',...
    [0.77 0.83 0.067 0.085],...
    'String',['~s^{',num2str(P52(1),'%.2f'),'\pm',num2str(E52(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H5,'textbox',...
    [0.77 0.36 0.067 0.085],...
    'String',['~s^{',num2str(P53(1),'%.2f'),'\pm',num2str(E53(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(H5,'textbox',...
    [0.33 0.36 0.067 0.085],...
    'String',['~s^{',num2str(P54(1),'%.2f'),'\pm',num2str(E54(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

ft = fit(log(mydist(:)),Y20(:),ftype1,'Weight',1./V20(:));
CI20 = confint(ft,0.68);
P20 = [ft.a,ft.b];
E20 = mean(abs(CI20-repmat(P20,2,1)));
p=plot(h20,x,exp(P20(2))*x.^P20(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(log(mydist(:)),Y20_global(:),ftype1,'Weight',1./V20_global(:));
CI20_global = confint(ft,0.68);
P20_global = [ft.a,ft.b];
E20_global = mean(abs(CI20_global-repmat(P20_global,2,1)));
p=plot(h20,x,exp(P20_global(2))*x.^P20_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
uistack(p,'bottom')

annotation(H20,'textbox',...
    [0.15 0.835 0.067 0.085],...
    'String',['~s^{',num2str(P20_global(1),'%.2f'),'\pm',num2str(E20_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% Compare

%P(Ooff) P(Poff) P(Pon) T(Poff) T(Pon) T(P) relaxation time
meanval = [57.8000000000000,0.523255813953488,0.175429726996966,0.301314459049545,5.51954022988505,12.9000000000000,23.4000000000000,2.07012993870273;...
    81.7000000000000,0.356373429084380,0.466487133453022,0.177139437462597,8.83209876543209,11.3000000000000,17.7000000000000,3.77385998539811;...
    87.5000000000000,0.650273224043715,0.233151183970856,0.116575591985428,7.50350877192982,13,15.1000000000000,3.58986779799169;...
    148.700000000000,0.834322829451463,0.106916324148779,0.0587608463997564,8.54672489082969,12.4000000000000,16,5.25409164776652;...
    189.600000000000,0.764133016627078,0.132185273159144,0.103681710213776,8.08387096774193,15.5000000000000,17.4000000000000,6.73814095177142;...
    595.100000000000,0.931463127834705,0.0477070384679993,0.0208298336972954,8.07916666666666,NaN,13.6000000000000,9.96348444050395];

stdval = [57.8000000000000,0.0607529725735618,0.0320528925062542,0.0513425727289546,1.11113472031996,2.94000000000000,5.55000000000000,0.249608414878848;...
    81.7000000000000,0.0446895619493190,0.0384603357886948,0.0352434491102140,1.16013708580682,2.44000000000000,2.70000000000000,0.530332900154374;...
    87.5000000000000,0.0464440234734158,0.0376845071723769,0.0313690539865890,1.35533123668975,3.76000000000000,3.04000000000000,0.346914528268444;...
    148.700000000000,0.0124924930534859,0.00846300616158664,0.00737875673382437,0.705287379927667,1.59000000000000,1.43000000000000,0.163538322016393;...
    189.600000000000,0.0299190477348723,0.0195513955706115,0.0221601995867767,1.19100438717683,3.25000000000000,2.73000000000000,0.333976511907366;...
    595.100000000000,0.0194376995828243,0.0145738984286442,0.0107645523676546,2.03984697070557,NaN,3.61000000000000,0.520591351687275];

mm = nanmean(meanval,1);
em = sqrt(nanvar(meanval,1,1)+nanmean(stdval.^2,1))./sqrt(sum(~isnan(meanval),1));

H6=figure(6);
set(H6,'position',[50 700 3*Wi 1*Le],'paperpositionmode','auto','color','w');
h61 = subplot(1,3,1,'parent',H6);
h62 = subplot(1,3,2,'parent',H6);
h63 = subplot(1,3,3,'parent',H6);
hold(h61,'on')
hold(h62,'on')
hold(h63,'on')
set(h61,'fontsize',ftsz,'linewidth',2,'tickdir','out')
set(h62,'fontsize',ftsz,'linewidth',2,'tickdir','out')
set(h63,'fontsize',ftsz,'linewidth',2,'tickdir','out')
xlabel(h61,'HMM P(O_{off})')
xlabel(h62,'HMM P(P_{off})')
xlabel(h63,'HMM P(P_{on})')
ylabel(h61,'data P(O_{off})')
ylabel(h62,'data P(P_{off})')
ylabel(h63,'data P(P_{on})')
xlim(h61,[0,1])
xlim(h62,[0,1])
xlim(h63,[0,1])
ylim(h61,[0,1])
ylim(h62,[0,1])
ylim(h63,[0,1])
box(h61,'on')
box(h62,'on')
box(h63,'on')

H7=figure(7);
set(H7,'position',[50 700 3*Wi 1*Le],'paperpositionmode','auto','color','w');
h71 = subplot(1,3,1,'parent',H7);
h72 = subplot(1,3,2,'parent',H7);
h73 = subplot(1,3,3,'parent',H7);
hold(h71,'on')
hold(h72,'on')
hold(h73,'on')
set(h71,'fontsize',ftsz,'linewidth',2,'tickdir','out')
set(h72,'fontsize',ftsz,'linewidth',2,'tickdir','out')
set(h73,'fontsize',ftsz,'linewidth',2,'tickdir','out')
xlabel(h71,'HMM T(P) (min)')
xlabel(h72,'HMM T(P_{off}) (min)')
xlabel(h73,'HMM T(P_{on}) (min)')
ylabel(h71,'data T(P) (min)')
ylabel(h72,'data T(P_{off}) (min)')
ylabel(h73,'data T(P_{on}) (min)')
xlim(h71,[0,50])
xlim(h72,[0,15])
xlim(h73,[0,40])
ylim(h71,[0,50])
ylim(h72,[0,15])
ylim(h73,[0,40])
box(h71,'on')
box(h72,'on')
box(h73,'on')

H8=figure(8);
set(H8,'position',[50 700 1.2*Wi 1.35*Le],'paperpositionmode','auto','color','w');
h81 = subplot(1,1,1,'parent',H8);
hold(h81,'on')
set(h81,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h81,'data \tau [min]')
ylabel(h81,'HMM T(O_{off}) [min]')
box(h81,'on')

H9=figure(9);
set(H9,'position',[50 700 1.2*Wi 1.35*Le],'paperpositionmode','auto','color','w');
h91 = subplot(1,1,1,'parent',H9);
hold(h91,'on')
set(h91,'fontsize',ftsz,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h91,'s (kb)')
ylabel(h91,'T(O_{off})/\tau')
box(h91,'on')
xlim(h91,[40,1000])
ylim(h91,[1,100])

plot(h61,[0,1],[0,1],'--k','linewidth',lnwi)
plot(h62,[0,1],[0,1],'--k','linewidth',lnwi)
plot(h63,[0,1],[0,1],'--k','linewidth',lnwi)

plot(h71,[0,50],[0,50],'--k','linewidth',lnwi)
plot(h72,[0,15],[0,15],'--k','linewidth',lnwi)
plot(h73,[0,40],[0,40],'--k','linewidth',lnwi)

X81 = nan(Nd,1);
Y91 = nan(Nd,1);
V91 = nan(Nd,1);
Y91_global = nan(Nd,1);
V91_global = nan(Nd,1);

M91 = nan(Nd,1);
S91 = nan(Nd,1);

for i=1:Nd
    errorbar(h61,exp(Y31_global(i)),meanval(i,2),stdval(i,2),stdval(i,2),EB31_global(i,1),EB31_global(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h62,exp(Y32_global(i)),meanval(i,3),stdval(i,3),stdval(i,3),EB32_global(i,1),EB32_global(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h63,exp(Y33_global(i)),meanval(i,4),stdval(i,4),stdval(i,4),EB33_global(i,1),EB33_global(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)

    errorbar(h61,exp(Y31(i)),meanval(i,2),stdval(i,2),stdval(i,2),EB31(i,1),EB31(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    errorbar(h62,exp(Y32(i)),meanval(i,3),stdval(i,3),stdval(i,3),EB32(i,1),EB32(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    errorbar(h63,exp(Y33(i)),meanval(i,4),stdval(i,4),stdval(i,4),EB33(i,1),EB33(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    
    errorbar(h71,exp(Y54(i)),meanval(i,7),stdval(i,7),stdval(i,7),EB54(i,1),EB54(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    errorbar(h72,exp(Y52(i)),meanval(i,5),stdval(i,5),stdval(i,5),EB52(i,1),EB52(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    errorbar(h73,exp(Y53(i)),meanval(i,6),stdval(i,6),stdval(i,6),EB53(i,1),EB53(i,2),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    
    errorbar(h81,meanval(i,8),exp(Y51_global(i)),EB51_global(i,1),EB51_global(i,2),stdval(i,8),stdval(i,8),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    errorbar(h81,meanval(i,8),exp(Y51(i)),EB51(i,1),EB51(i,2),stdval(i,8),stdval(i,8),'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    
    ey = stdval(i,8);
    ex = mean([EB51_global(i,1),EB51_global(i,2)]);
    y = (meanval(i,8));
    x = exp(Y51_global(i));
    ee = sqrt(ex^2/y^2+ey^2*x^2/y^4);
    
    M91(i) = x/y;
    S91(i) = ee;
    
    errorbar(h91,mydist(i),x/y,ee,'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
    Y91_global(i) = log(x/y);
    V91_global(i) = (0.5*(log(ee+x/y)-log(x/y-ee)))^2;
    
    ex = mean([EB51(i,1),EB51(i,2)]);
    x = exp(Y51(i));
    ee = sqrt(ex^2/y^2+ey^2*x^2/y^4);
    
    errorbar(h91,mydist(i),x/y,ee,'o','linewidth',lnwi,'markersize',mksz,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'capsize',0)
    Y91(i) = log(x/y);
    V91(i) = (0.5*(log(ee+x/y)-log(x/y-ee)))^2;
    
    X81(i) = log(meanval(i,8));
end

errorbar(h71,Tm_global(4),mm(7),em(7),em(7),Tm_global(4)-Tl_global(4),Tu_global(4)-Tm_global(4),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
errorbar(h72,Tm_global(2),mm(5),em(5),em(5),Tm_global(2)-Tl_global(2),Tu_global(2)-Tm_global(2),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)
errorbar(h73,Tm_global(3),mm(6),em(6),em(6),Tm_global(3)-Tl_global(3),Tu_global(3)-Tm_global(3),'o','linewidth',lnwi,'markersize',mksz,'color',[0,0,0],'markerfacecolor',[0,0,0],'capsize',0)

mr = nanmean(M91,1);
er = sqrt(nanvar(M91,1,1)+nanmean(S91.^2,1))./sqrt(sum(~isnan(M91),1));

p=plot(h91,[40,1000],mr*[1,1],'-k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h91,[40,1000],(mr+er)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')
p=plot(h91,[40,1000],(mr-er)*[1,1],':k','linewidth',lnwi);
uistack(p,'bottom')

x8 = linspace(1,20,100);

ft = fit(X81(:),Y51(:),ftype1,'Weight',1./V51(:));
CI81 = confint(ft,0.68);
P81 = [ft.a,ft.b];
E81 = mean(abs(CI81-repmat(P81,2,1)));
p=plot(h81,x8,exp(P81(2))*x8.^P81(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
uistack(p,'bottom')

ft = fit(X81(:),Y51_global(:),ftype1,'Weight',1./V51_global(:));
CI81_global = confint(ft,0.68);
P81_global = [ft.a,ft.b];
E81_global = mean(abs(CI81_global-repmat(P81_global,2,1)));
p=plot(h81,x8,exp(P81_global(2))*x8.^P81_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
uistack(p,'bottom')

annotation(H8,'textbox',...
    [0.22 0.82 0.067 0.085],...
    'String',['~\tau^{',num2str(P81_global(1),'%.2f'),'\pm',num2str(E81_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

x9 = logspace(log10(3e1),log10(1e3),100);

ft = fit(log(mydist(:)),Y91(:),ftype1,'Weight',1./V91(:));
CI91 = confint(ft,0.68);
P91 = [ft.a,ft.b];
E91 = mean(abs(CI91-repmat(P91,2,1)));
%p=plot(h91,x9,exp(P91(2))*x9.^P91(1),'-','linewidth',lnwi,'color',[0.5,0.5,0.5]);
%uistack(p,'bottom')

ft = fit(log(mydist(:)),Y91_global(:),ftype1,'Weight',1./V91_global(:));
CI91_global = confint(ft,0.68);
P91_global = [ft.a,ft.b];
E91_global = mean(abs(CI91_global-repmat(P91_global,2,1)));
%p=plot(h91,x9,exp(P91_global(2))*x9.^P91_global(1),'-','linewidth',lnwi,'color',[0,0,0]);
%uistack(p,'bottom')

annotation(H9,'textbox',...
    [0.61 0.82 0.067 0.085],...
    'String',['~s^{',num2str(P91_global(1),'%.2f'),'\pm',num2str(E91_global(1),'%.2f'),'}'],...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

