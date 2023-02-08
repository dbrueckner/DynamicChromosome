function [ F, Ln, L ] = Forward( Po, Pt, P0, tt )
%Forward algorithm
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

Ns = size(P0,1);
Nd = size(Po,2);
Ln = zeros(Nd,1);

F  = zeros(Ns,Nd);
Fi = Po(:,1).*P0;

Ln(1) = sum(Fi);
Fi = Fi / Ln(1);
Ln(1) = log(Ln(1));

F(:,1) = Fi;

for i=2:Nd
    v = tt(i-1);
    Fj = Po(:,i).* (Pt{v} * Fi);
    Ln(i) = sum(Fj);
    Fj = Fj / Ln(i);
    Ln(i) = log(Ln(i));
    
    F(:,i) = Fj;
    Fi = Fj;
end

L = sum(Ln);

end

