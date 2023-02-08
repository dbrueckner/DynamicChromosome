%MainInferenceHomie
%   Perform HMM inference on homie data
%
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

clc
clear
close all

datapath = '/Users/bzoller/PostDoc/Writing/DavidLeo_Paper/Processing/LinesData/';

timeRes = 0.46;

if ~exist('FullDataSet.mat','file')
    %%% Make input Data structure
    s = [58,82,88,149,190,595];
    blue = cell(6,2);
    green = cell(6,2);
    bgVec = cell(6,2);
    bgDist = cell(6,2);
    red = cell(6,2);
    
    for i=1:6
        [blue_i,green_i,bgVec_i,bgDist_i,red_i] = ReadDataHomie([datapath,'data_line',num2str(i-1),'.csv']);
        
        l=length(blue_i);
        blue(i,1:l) = blue_i;
        green(i,1:l) = green_i;
        bgVec(i,1:l) = bgVec_i;
        bgDist(i,1:l) = bgDist_i;
        red(i,1:l) = red_i;
    end
    save('FullDataSet.mat','s','blue','green','bgVec','bgDist','red')
else
    %%% Load input Data structure
    load('FullDataSet.mat','s','blue','green','bgVec','bgDist','red')
end

%MCMCmode = 'hierachichal';
MCMCmode = 'single';

switch MCMCmode
    case 'hierachichal'
        % all genomic separation together
        vs = 1:6;
    case 'single'
        % single genomic separation
        vs = 4;
end

%% MCMC Sampling
count = 50000;
MCRepeats = 1;
Nblock = length(vs);

if Nblock>1
    %PoffD=PonD=PD
    %global:b1,f2,b2,b3,PD,red1,red2,red3
    MCMCStruct.var.global = 8; 
    %f1,OpenD
    MCMCStruct.var.local = 2;
else
    %fixed:b3,red1,red2,red3
    Xfixed = [0,0.114619368133244,0.160319751888197,0.951851816181802];
    %total:b1,f2,b2,PoffD=PonD,f1,OpenD
    MCMCStruct.var.global = 6;
    MCMCStruct.var.local = 0;
end

%%% input:f1,b1,f2,b2,(b3),OpenD,PoffD,PonD,red1,red2,red3
%%% initial condition

%global:b1,f2,b2,b3,PoffD,PonD,red1,red2,red3
%PoffD=PonD=PD
x0global = [rand(), rand(), rand(), rand(), 400, .1, .1, 1];

%local:f1,OpenD
x0local = [rand(1,6);[600, 660, 720, 780, 840, 900]];
x0 = x0local(:,vs);
x0 = [x0global, x0(:)'];

if Nblock==1
    x0 = [rand(), rand(), rand(), 400, rand(), x0local(2,vs)];
end

%%% configure likelihood and prior
if Nblock>1
    %global:b1,f2,b2,(b3),PoffD,PonD,red1,red2,red3
    %local:f1,OpenD
    for i=1:length(vs)
        %PoffD=PonD=PD
        MCMCStruct.block(i).loglikelihood = @(X,Y) LogLikelihood(bgDist(vs(i),:),red(vs(i),:),[Y(1),X(1:4)],[Y(2),X(5),X(5)],X(6:8),timeRes);
        MCMCStruct.block(i).logprior = @(X,Y) LogPrior(Y(1),Y(2),nan);
    end    
    %PoffD=PonD=PD
    MCMCStruct.global.logprior = @(X) LogPrior(X(1:4),X(5),X(6:8));
else
    %Xfixed
    MCMCStruct.block(1).loglikelihood = @(X,Y) LogLikelihood(bgDist(vs(1),:),red(vs(1),:),[X(5),X(1:3),Xfixed(1)],[X(6),X(4),X(4)],Xfixed(2:4),timeRes);
    MCMCStruct.block(1).logprior = @(X,Y) 0;
    MCMCStruct.global.logprior = @(X) LogPrior([X(5),X(1:3)],[X(4),X(6)],nan);
end

%%% MCMC run
X = nan(MCRepeats*count,length(x0));
L = nan(MCRepeats*count,Nblock+1);
A = nan(MCRepeats*count,Nblock+1);
H = nan(MCRepeats*count,1);

tic
for mcCount = 1:MCRepeats
    [ X((1:count)+(mcCount-1)*count,:),L((1:count)+(mcCount-1)*count,:),A((1:count)+(mcCount-1)*count,:),H((1:count)+(mcCount-1)*count) ] = MCMC_block_sampler( MCMCStruct, x0, [], count, 0.234);
    x0 = X(mcCount*count,:);
end
toc

%%% save MCMC run
switch MCMCmode
    case 'hierachichal'
        save('Hierachichal_MCMC.mat');
    case 'single'
        save(['Single',num2str(s(vs)),'_MCMC.mat']);
end

%% Posterior decoding
clc
close all
Wi = 500;
Le = 375;

plotTrajectories = false;
burnin = 2e4;

%%% load MCMC run
load('Hierachichal_MCMC.mat');
X = [X(:,1:5),X(:,5),X(:,6:end)];
Xest = mean(X(burnin:end,:));

state = {};
for i = 1:length(vs)
    %empty cells checked
    trackLength = sum(cellfun(@(x) ~isempty(x),bgDist(vs(i),:)));
    for j = 1:trackLength
        currentTrackDist = bgDist{vs(i),j};
        currentTrackRed = red{vs(i),j};
        %%% input
        % f1,b1,f2,b2
        % OpenD,PoffD,PonD
        % red1,red2,red3
        s = PosteriorDecoding({currentTrackDist},{currentTrackRed},[Xest(10+2*(i-1)),Xest(1:4)],[Xest(11+2*(i-1)),Xest(5:6)],Xest(7:9),timeRes);
        v = ViterbiDecoding({currentTrackDist},{currentTrackRed},[Xest(10+2*(i-1)),Xest(1:4)],[Xest(11+2*(i-1)),Xest(5:6)],Xest(7:9),timeRes);
        
        [~,smax] = max(s(:,3:5),[],2);
        smax(isnan(s(:,1))) = nan;
        state{i,j} = smax-1;
        
        if(plotTrajectories)
            nt = length(bgDist{i,j});
            tt = (0:(nt-1))*timeRes;
            
            Hi=figure(134);
            set(Hi,'position',[50 700 1.5*Wi 2*Le],'paperpositionmode','auto','color','w');
            
            hi1=subplot(311,'parent',Hi);
            hold(hi1,'on')
            box(hi1,'on')
            plot(hi1,tt,currentTrackDist,'-','linewidth',3)
            set(hi1,'fontsize',24,'linewidth',2,'tickdir','out')
            ylim(hi1,[0,2000]);
            title(hi1,'blue-green distances')
            ylabel(hi1,'R (nm)')
            xlim(hi1,[tt(1),tt(end)]);
            set(hi1,'xtick',0:10:40)
            
            hi2=subplot(312,'parent',Hi);
            hold(hi2,'on')
            box(hi2,'on')
            plot(hi2,tt,currentTrackRed,'-r','linewidth',3)
            set(hi2,'fontsize',24,'linewidth',2,'tickdir','out')
            ylim(hi2,[0,4]);
            title(hi2,'red channel intensity')
            ylabel(hi2,'I (a.u.)');
            xlim(hi2,[tt(1),tt(end)])
            set(hi2,'xtick',0:10:40)
            
            hi3=subplot(313,'parent',Hi);
            hold(hi3,'on')
            box(hi3,'on')
            plot(hi3,tt,smax,'-m','linewidth',3)
            plot(hi3,tt,v,'--g','linewidth',3)
            set(hi3,'fontsize',24,'linewidth',2,'tickdir','out','ytick',1:3,'yticklabel',{'O_{off}','P_{off}','P_{on}'})
            title(hi3,'state')
            ylim(hi3,[0.5,3.5]);
            xlim(hi3,[tt(1),tt(end)])
            set(hi3,'xtick',0:10:40,'ytick',1:3)
            
            pause()
            close(Hi)
        end
    end
end

results = X;
save('FullDataSet','s','blue','green','bgVec','bgDist','red','state','results')

%% Likelihood and prior

function [ L ] = LogLikelihood( traces_d, traces_r, kinetics, distances, intensities, dt )

f1 = kinetics(1);
b1 = kinetics(2);
f2 = kinetics(3);
b2 = kinetics(4);
b3 = kinetics(5);

M = [-f1, b1, b3;
    f1, -b1-f2, b2;
    0, f2, -b2-b3];
Ns = size(M,1);

%1,2,3
P0 = SteadyState( M );
%safety check
P0 = abs(P0);
P0 = P0/sum(P0);

%50 missing time points in a row should be enough
Pt = cell(50,1);
for j = 1:50
    Pt{j} = expm(M*(dt*j));
end

emptyTraces = cellfun(@(X) isempty(X),traces_d);
traces_d_I = traces_d(~emptyTraces);
traces_r_I = traces_r(~emptyTraces);

LL = nan(length(traces_d_I),1);
for j = 1:length(traces_d_I)
    D = traces_d_I{j};
    R = traces_r_I{j};
    
    %let's remove bad values beforehand
    mask = ~(D==0 | R==0 | isnan(D) | isnan(R));
    counter = 1:length(D);
    tt = diff(counter(mask));
    
    Po = ObservationMatrix( D(mask), R(mask), distances, intensities, Ns );
    [~,~,LL(j)] = Forward( Po, Pt, P0, tt );
end
L = nansum(LL);

end

function out = LogPrior(kinetics,distances,intensities)

kMask = kinetics > 5 | kinetics <= 1e-6;
kinetics(kMask) = -1e7;
kinetics(~kMask) = -log(kinetics(~kMask));

dMask = distances > 2500 | distances <= 200;
distances(dMask) = -1e7;
distances(~dMask) = -log(distances(~dMask));

rMask = intensities > 5 | intensities <= 1e-6;
intensities(rMask) = -1e7;
intensities(~rMask) = -log(intensities(~rMask));

out = nansum(kinetics)+nansum(distances)+nansum(intensities);

end
