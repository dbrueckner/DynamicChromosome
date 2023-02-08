function [ ss_state ] = SteadyState( M )
%SteadyState
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

M(end,:) = ones(1,size(M,1));
b = zeros(size(M,1),1);
b(end) = 1;
ss_state = M \ b;

end

