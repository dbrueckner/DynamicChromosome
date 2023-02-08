function S = Viterbi( Po, Pt, P0, tt )
%Viterbi algorithm
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

Ns = size(P0,1);
Nd = size(Po,2);

Ki = zeros(Ns,Nd);
Vi = Po(:,1).*P0;

for i=2:Nd
    v = tt(i-1);
    [Mi,Ki(:,i)] = max(Pt{v} .* repmat(Vi',Ns,1),[],2);  
    Vi = Po(:,i).*Mi;
    Vi = Vi/sum(Vi);
end

[~,k] = max(Vi);

S = zeros(1,Nd);
S(end) = k;

for i=Nd:-1:2
    S(i-1) = Ki(S(i),i);
end

end

