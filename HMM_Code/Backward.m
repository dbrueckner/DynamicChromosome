function [ B, Ln, L ] = Backward( Po, Pt, P0, tt )
%Backward algorithm
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

Ns = size(P0,1);
Nd = size(Po,2);
Ln = zeros(Nd,1);

B  = zeros(Ns,Nd);
Bi = ones(Ns,1)/Ns;
Ln(end) = log(Ns);

B(:,end) = Bi;

for i=(Nd-1):-1:1
    v = tt(i);
    Bj = (Pt{v})' * (Po(:,i+1).* Bi);
    Ln(i) = sum(Bj);
    Bj = Bj / Ln(i);
    Ln(i) = log(Ln(i));
    
    B(:,i) = Bj;
    Bi = Bj;
end

Bj = P0 .* Po(:,1) .* Bi;
L = log(sum(Bj)) + sum(Ln);

end

