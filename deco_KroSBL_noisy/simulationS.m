function [resultS] = simulationS(y,A1,A2,A3,A,N,R_max,x,SNR,epi)

% SNR noise
signal_power = (norm(y)^2)/length(y);
noise_var = (signal_power)/(10^(SNR/10));
noise = sqrt(noise_var)*randn(size(y));
y = y + noise;

% implementation of each algorithm
%% classic SBL
[metrics_csbl] = classicSBL(y,A,N,R_max,x);
%% AM_KroSBL
[metrics_am] = am_kroSBL(y,A1,A2,A3,A,N,R_max,x);
%% SVD_KroSBL
[metrics_svd] = svd_kroSBL(y,A1,A2,A3,A,N,R_max,x);
%% deco_KroSBL
[metrics_deco] = deco_kroSBL(y,A1,A2,A3,N,R_max,x);
%% OMP
[metrics_omp] = omp2(y,A,x,epi);
% %% fista
% [metrics_fista] = deco_fista_lasso(y,A1,A2,A3,N,x);
%% deco_omp
[metrics_deco_omp] = deco_OMP(y,A1,A2,A3,N,x,epi);
%% data saving
result = cell(7,1);
result{1} = metrics_csbl;
result{2} = metrics_am;
result{3} = metrics_svd;
result{4} = metrics_deco;
result{5} = metrics_omp;
result{6} = metrics_deco_omp;
result{7} = noise_var;
resultS = result;