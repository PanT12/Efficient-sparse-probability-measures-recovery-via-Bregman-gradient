function [number,mean,deviance_mean,correlation] = importfile(filename, dataLines)
% The format of these data files is:
% number of assets (N)
% for each asset i (i=1,...,N):
%    mean return, standard deviation of return
% for all possible pairs of assets:
%    i, j, correlation between asset i and asset j

% the following codes are generated in Matlab.
%% setting
if nargin < 2
    dataLines = [1, 528];
end

opts = delimitedTextImportOptions("NumVariables", 4);

opts.DataLines = dataLines;
opts.Delimiter = " ";

opts.VariableNames = ["VarName1", "VarName2", "VarName3", "Var4"];
opts.SelectedVariableNames = ["VarName1", "VarName2", "VarName3"];
opts.VariableTypes = ["double", "double", "double", "string"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

opts = setvaropts(opts, "Var4", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var4", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "VarName1", "TrimNonNumeric", true);
opts = setvaropts(opts, "VarName1", "ThousandsSeparator", ",");

port1 = readtable(filename, opts);


port_data = table2array(port1);
number = port_data(1,1);
mean = port_data(2:(number+1),1);
deviance_mean = port_data(2:(number+1),2);
% construct correlation matrix
port_data = port_data((number+2):end,3);
correlation = ones(number,number);
a = 1;
for i = 1:number
    correlation(i,i:end) = port_data(a:(a+number-i));
    a = a+number-i + 1;
end
correlation = triu(correlation,0) + tril(correlation',-1);
end