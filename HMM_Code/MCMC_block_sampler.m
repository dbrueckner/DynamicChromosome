function [ X,L,A,H ] = MCMC_block_sampler( Likeli_Struct, Xprev, Lprev, N_iter, ar, glob_prop, local_prop )
%MCMC_block_sampler
%   Adapative block MCMC sampler
%
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

if nargin < 7
    glob_prop  = 'lognormal';
    local_prop = 'lognormal';
end

Nvar_global = Likeli_Struct.var.global;
Nvar_local = Likeli_Struct.var.local;
Nblock = length(Likeli_Struct.block);
% loglikelihood = Likeli_Struct.block(i).loglikelihood;
% function handle @(X,Y) where X =
% [x(1),...,x(Nvar_global)],Y = [y(1),...,y(Nvar_local)]
% logprior = Likeli_Struct.block(i).logprior;
% function handle @(X,Y) here X =
% [x(1),...,x(Nvar_global)],Y = [y(1),...,y(Nvar_local)]
% logprior = Likeli_Struct.global.logprior;
% function handle @(X) where X =
% [x(1),...,x(Nvar_global)]

Nvar_tot = Nvar_global + Nblock*Nvar_local;
X = zeros(N_iter,Nvar_tot);
L = zeros(N_iter,Nblock+1);
A = zeros(N_iter,Nblock+1);
H = zeros(N_iter,1);

if isempty(Lprev)
    Lprev = zeros(Nblock+1,1);
    Xprev_global = Xprev(1:Nvar_global);
    
    for i=1:Nblock
        Xprev_local = Xprev(Nvar_global+(i-1)*Nvar_local+(1:Nvar_local));
        
        loglikelihood = Likeli_Struct.block(i).loglikelihood;
        logprior = Likeli_Struct.block(i).logprior; % only for local variables
        
        Lprev(i+1) = loglikelihood(Xprev_global,Xprev_local) + logprior(Xprev_global,Xprev_local);
    end
    
    Lprev(1) = sum(Lprev(2:end)) + Likeli_Struct.global.logprior(Xprev_global);
end

move_set = {'global','local'};
if Nvar_local == 0
    pglobal = 1.01;
else
    pglobal = 1/3;
end

% For Adaptation
l = cell(Nblock+1,1);
S = cell(Nblock+1,1);
M = cell(Nblock+1,1);

switch glob_prop
    case 'lognormal'
        M{1} = log(Xprev(1:Nvar_global));
    case 'normal'
        M{1} = Xprev(1:Nvar_global);
end
l{1} = 1;
S{1} = eye(Nvar_global)/10;

for i=1:Nblock
    Xprev_local = Xprev(Nvar_global+(i-1)*Nvar_local+(1:Nvar_local));
    switch local_prop
        case 'lognormal'
            M{1+i} = log(Xprev_local);
        case 'normal'
            M{1+i} = Xprev_local;
    end
    l{1+i} = 1;
    S{1+i} = eye(Nvar_local)/10;
end

for i=1:N_iter
    r = rand(1);
    if r < pglobal
        move = move_set{1};
    else
        move = move_set{2};
    end
    
    switch move
        case 'global'
            [~,err] = cholcov(l{1}*S{1});
            if err ~= 0
                fprintf('error: global proposal covariance diverged\n');
                return
            end
            
            H(i) = 1;
            n = sum(H(1:i));
            [Xprop,Lprop,Accept,l,S,M] = Acceptance_global(Likeli_Struct,Xprev,Lprev,l,S,M,n,ar,glob_prop);
            
        case 'local'
            H(i) = 0;
            n = sum(~H(1:i));
            [Xprop,Lprop,Accept,l,S,M] = Acceptance_local(Likeli_Struct,Xprev,Lprev,l,S,M,n,ar,local_prop);
            
    end
    
    X(i,:) = Xprop;
    L(i,:) = Lprop;
    A(i,:) = Accept;
    
    if mod(i,1)==0
        fprintf('iter %i: %f %f \n',i,Lprop(1),Lprev(1));
    end
    
    Xprev = Xprop;
    Lprev = Lprop; 
end

end

function [Xprop,Lprop,Accept,l,S,M] = Acceptance_global(Likeli_Struct,Xprev,Lprev,l,S,M,n,ar,glob_prop)
    
    Nvar_global = Likeli_Struct.var.global;
    Nvar_local = Likeli_Struct.var.local;
    Nblock = length(Likeli_Struct.block);

    Xprop = zeros(size(Xprev));
    Lprop = zeros(size(Lprev));
    Accept = zeros(size(Lprev));
    
    Xprop((Nvar_global+1):end) = Xprev((Nvar_global+1):end);
    Xprev_global = Xprev(1:Nvar_global);
    switch glob_prop
        case 'lognormal'
            Xprop_global = exp(mvnrnd(log(Xprev_global),l{1}*S{1},1));
        case 'normal'
            Xprop_global = mvnrnd(log(Xprev_global),l{1}*S{1},1);
    end
    
    Lprev_global = Lprev(1);
    Lprop_local = zeros(Nblock,1);
    
    for i=1:Nblock
        Xprev_local = Xprev(Nvar_global+(i-1)*Nvar_local+(1:Nvar_local));
        
        loglikelihood = Likeli_Struct.block(i).loglikelihood;
        logprior = Likeli_Struct.block(i).logprior; % only for local variables
        
        Lprop_local(i) = loglikelihood(Xprop_global,Xprev_local) + logprior(Xprop_global,Xprev_local);
    end
    
    Lprop_global = sum(Lprop_local) + Likeli_Struct.global.logprior(Xprop_global);
    
    switch glob_prop
        case 'lognormal'
            Qprop_global = -sum(log(Xprop_global));
            Qprev_global = -sum(log(Xprev_global));
        case 'normal'
            Qprop_global = 0;
            Qprev_global = 0;
    end
    
    aprob = min(1,exp(Lprop_global+Qprev_global-Lprev_global-Qprop_global));
    
    r = rand(1);
    
    if r < aprob
        Xprop(1:Nvar_global) = Xprop_global;
        Lprop(1) = Lprop_global;
        Lprop(2:end) = Lprop_local;
        Accept(1) = 1;
    else
        Xprop_global = Xprev_global;
        Xprop(1:Nvar_global) = Xprop_global;
        Lprop(1) = Lprev_global;
        Lprop(2:end) = Lprev(2:end);
        Accept(1) = 0;
    end
    
    % Adaptation
    y = 1/(3*sqrt(n));
    l{1} = exp(log(l{1}) + y*(aprob - ar));
    switch glob_prop
        case 'lognormal'
            S{1} = S{1} + y*((log(Xprop_global) - M{1})'*(log(Xprop_global) - M{1}) - S{1});
            M{1} = M{1} + y*(log(Xprop_global) - M{1});
        case 'normal'
            S{1} = S{1} + y*((Xprop_global - M{1})'*(Xprop_global - M{1}) - S{1});
            M{1} = M{1} + y*(Xprop_global - M{1});
    end

end

function [Xprop,Lprop,Accept,l,S,M] = Acceptance_local(Likeli_Struct,Xprev,Lprev,l,S,M,n,ar,local_prop)

    Nvar_global = Likeli_Struct.var.global;
    Nvar_local = Likeli_Struct.var.local;
    Nblock = length(Likeli_Struct.block);
    
    Xprop = zeros(size(Xprev));
    Lprop = zeros(size(Lprev));
    Accept = zeros(size(Lprev));
    
    Xprev_global = Xprev(1:Nvar_global);
    Xprop(1:Nvar_global) = Xprev_global;
    
    for i=1:Nblock
        loglikelihood = Likeli_Struct.block(i).loglikelihood;
        logprior = Likeli_Struct.block(i).logprior; % only for local variables
        
        Xprev_local = Xprev(Nvar_global+(i-1)*Nvar_local+(1:Nvar_local));
        Lprev_local = Lprev(1+i);
        
        switch local_prop
            case 'lognormal'
                Xprop_local = exp(mvnrnd(log(Xprev_local),l{i+1}*S{i+1},1));
            case 'normal'
                Xprop_local = mvnrnd(log(Xprev_local),l{i+1}*S{i+1},1);
        end
        Lprop_local = loglikelihood(Xprev_global,Xprop_local) + logprior(Xprev_global,Xprop_local);
        
        switch local_prop
            case 'lognormal'
                Qprop_local = -sum(log(Xprop_local));
                Qprev_local = -sum(log(Xprev_local));
            case 'normal'
                Qprop_local = 0;
                Qprev_local = 0;
        end
        
        aprob = min(1,exp(Lprop_local+Qprev_local-Lprev_local-Qprop_local));
        
        r = rand(1);
        
        if r < aprob
            Xprop(Nvar_global+(i-1)*Nvar_local+(1:Nvar_local)) = Xprop_local;
            Lprop(1+i) = Lprop_local;
            Accept(1+i) = 1;
        else
            Xprop_local = Xprev_local;
            Xprop(Nvar_global+(i-1)*Nvar_local+(1:Nvar_local)) = Xprop_local;
            Lprop(1+i) = Lprev_local;
            Accept(1+i) = 0;
        end
        
        % Adaptation
        y = 1/(3*sqrt(n));
        l{i+1} = exp(log(l{i+1}) + y*(aprob - ar));
        switch local_prop
            case 'lognormal'
                S{i+1} = S{i+1} + y*((log(Xprop_local) - M{i+1})'*(log(Xprop_local) - M{i+1}) - S{i+1});
                M{i+1} = M{i+1} + y*(log(Xprop_local) - M{i+1});                
            case 'normal'
                S{i+1} = S{i+1} + y*((Xprop_local - M{i+1})'*(Xprop_local - M{i+1}) - S{i+1});
                M{i+1} = M{i+1} + y*(Xprop_local - M{i+1});
        end
        
    end
       
    Lprop(1) = sum(Lprop(2:end)) + Likeli_Struct.global.logprior(Xprev_global);

end

