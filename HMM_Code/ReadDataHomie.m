function [blue,green,bgVec,bgDist,red] = ReadDataHomie(filename)
%ReadDataHomie
%   Copyright (c) 2023, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

%Setup the Import Options and import the bgDist and red intensities
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["cell_id", "time_point", "x_blue", "y_blue", "z_blue", "x_green", "y_green", "z_green", "x_Rij", "y_Rij", "z_Rij", "red", "Var13"];
opts.SelectedVariableNames = ["cell_id", "time_point", "x_blue", "y_blue", "z_blue", "x_green", "y_green", "z_green", "x_Rij", "y_Rij", "z_Rij", "red"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var13", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var13", "EmptyFieldRule", "auto");

% Import the bgDist and red intensities
mydata = readtable(filename, opts);

mydata = mydata{:,:};
Id = mydata(:,1)+1;

Nd = max(Id);
blue = cell(1,2);
green = cell(1,2);
bgVec = cell(1,2);
bgDist = cell(1,2);
red = cell(1,2);

for i=1:Nd
    Ii = Id == i;
    blue{1,i} = mydata(Ii,3:5);
    green{1,i} = mydata(Ii,6:8);
    bgvec = mydata(Ii,9:11);
    bgVec{1,i} = bgvec;
    bgDist{1,i} = sqrt(sum(bgvec.^2,2));
    red{1,i} = mydata(Ii,12);
end

end

