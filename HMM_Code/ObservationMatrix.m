function [ Po ] = ObservationMatrix( D, R, distances, intensities, Ns )
%ObservationMatrix
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

pobs = @(r,na) (3/(2*pi*na^2))^(3/2)*4*pi*r.^2.*exp(-3*r.^2/(2*na^2));
pobs_R = @(x,b) x./b^2 .*exp(-x.^2./(2.*b^2));

Nd = numel(D);
Po = zeros(Ns,Nd);

for i=1:Nd
    for j=1:Ns
        Po(j,i) = pobs(D(i),distances(j))*pobs_R(R(i),intensities(j));
    end
end
end
