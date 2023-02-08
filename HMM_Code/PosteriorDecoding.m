function [ PS ] = PosteriorDecoding( traces_d, traces_r, kinetics, distances, intensities, dt )
%PosteriorDecoding
%   Perform posterior decoding of time traces
%
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

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

Nt = length(traces_d);
PS = cell(Nt,1);

for j=1:Nt
    D = traces_d{j};
    R = traces_r{j};
    
    Nd = numel(D);
    PSj = nan(Nd,2+Ns);
    
    %let's remove bad values beforehand
    mask = ~(D==0 | R==0 | isnan(D) | isnan(R));
    counter = 1:length(D);
    tt = diff(counter(mask));
    
    Po = ObservationMatrix( D(mask), R(mask), distances, intensities, Ns );
    [ F,Lf_norm,Lf ] = Forward( Po, Pt, P0, tt );
    [ B,Lb_norm,Lb ] = Backward( Po, Pt, P0, tt );    
    
    flog_norm = cumsum(Lf_norm);
    blog_norm = cumsum(Lb_norm,'reverse');
    
    % Lf and Lb should give the same value
    [Lf,Lb]

    L = Lf;
    
    ei = 0:(Ns-1);
    ei = ei';
    
    Ni = sum(mask);
    v = 1:Nd;
    v = v(mask);
    for i=1:Ni
        Ps = F(:,i).*B(:,i) * exp(flog_norm(i) + blog_norm(i) - L);
        
        m1_e = sum(ei.*Ps);
        m2_e = sum(ei.^2.*Ps);
        s2_e = m2_e - m1_e^2;
        s1_e = sqrt(s2_e);
        
        PSj(v(i),:) = [m1_e,s1_e,Ps(:)'];
    end
    
    PS{j} = PSj;
end

if Nt==1
    PS = PSj;
end

end

