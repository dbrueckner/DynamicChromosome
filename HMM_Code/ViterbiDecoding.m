function VS = ViterbiDecoding( traces_d, traces_r, kinetics, distances, intensities, dt )
%ViterbiDecoding
%   Perform viterbi decoding of time traces
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
VS = cell(Nt,1);

for j=1:Nt
    D = traces_d{j};
    R = traces_r{j};
    
    Nd = numel(D);
    VSj = nan(Nd,1);
    
    %let's remove bad values beforehand
    mask = ~(D==0 | R==0 | isnan(D) | isnan(R));
    counter = 1:length(D);
    tt = diff(counter(mask));
    
    Po = ObservationMatrix( D(mask), R(mask), distances, intensities, Ns ); 
    S  = Viterbi( Po, Pt, P0, tt );
    
    ss = 1:Ns;
    
    Ni = sum(mask);
    v = 1:Nd;
    v = v(mask);
    for i=1:Ni
        VSj(v(i)) = ss(S(i));
    end
    
    VS{j} = VSj;
end

if Nt==1
    VS = VSj;
end

end

