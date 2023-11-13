% merge different files
clear
clc
x = load('compare_11.mat');
y = load('compare_37.mat');

vrs = fieldnames(x);
% Concatenate data

for i = 1:8
 x.('error'){i,2} = x.('error'){i,2} + y.('error'){i,2};
 x.('time'){i,2} = x.('time'){i,2} + y.('time'){i,2};
end

save('resultsss','-struct','x')