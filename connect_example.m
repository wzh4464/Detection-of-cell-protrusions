%% 
clear;
load('ssets.mat');
[disjoint_set,size_DS]=union_zero_fast(best,neib);
