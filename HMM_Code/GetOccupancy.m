function P = GetOccupancy( kinetics )
%GetOccupancy
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

Ns = size(kinetics,1);
P = zeros(Ns,3);

for i=1:Ns
    f1 = kinetics(i,1);
    b1 = kinetics(i,2);
    f2 = kinetics(i,3);
    b2 = kinetics(i,4);
    b3 = kinetics(i,5);
    
    M = [-f1, b1, b3;
        f1, -b1-f2, b2;
        0, f2, -b2-b3];
    
    %1,2,3
    P(i,:) = SteadyState( M );
end

end

