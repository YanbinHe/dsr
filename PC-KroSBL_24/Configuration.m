% fix random seed
rng(42)
%% simluation on the SBL based IRS channel estimation with angle spread
% predefined values
ms_ante = 6; % the number of antennas at user equippment side
bs_ante = 16; % the number of antennas at bs
irs_ele = 16^2; % the number of irs elements
Res1 = 18; % the number of atoms in the dicitionay
timeslots = 1e6;
modulating_scheme = 8;
epi = 0.1;

K1 = 4; % the number of overhead w.r.t pilot signals
K2 = 10; % the number of overhead w.r.t reflecting pattern

angle_spread = 3; % the spread of angles in the domain of cos
angle = linspace(-1,1-2/Res1,Res1); 

SNRl = [5,10,15,20,25,30]; % in dB
SNRl_10 = 10.^(SNRl/10);

numItr = 150; % the number of alternative optimization iterations
klevel = 50;

symbolmatrix = repmat((qammod((0:modulating_scheme-1)', modulating_scheme)),1,timeslots);
% generate pilot signals and irs pattern
X = 1/sqrt(ms_ante)*pilot_gen(ms_ante,max(K1));
irs_pattern = 1/sqrt(irs_ele)*((rand(irs_ele,max(K2)) > 0.5)*2-1); 
%% simulation settings
II = 1;
JJ = 6;
AVG = 50;% the number of trials
%%
% results
% nmse
error = cell(5,2);
error{1,1} = 'SVD-KroSBL';
error{1,2} = zeros(JJ,AVG);
error{2,1} = 'AM-KroSBL';
error{2,2} = zeros(JJ,AVG);
error{3,1} = 'dKroSBL';
error{3,2} = zeros(JJ,AVG);
error{4,1} = 'classicSBL';
error{4,2} = zeros(JJ,AVG);
error{5,1} = 'OMP';
error{5,2} = zeros(JJ,AVG);
error{6,1} = 'dOMP';
error{6,2} = zeros(JJ,AVG);
error{7,1} = 'dFISTA';
error{7,2} = zeros(JJ,AVG);


% support recovery
supprecovery = cell(5,2);
supprecovery{1,1} = 'SVD-KroSBL';
supprecovery{1,2} = zeros(JJ,AVG);
supprecovery{2,1} = 'AM-KroSBL';
supprecovery{2,2} = zeros(JJ,AVG);
supprecovery{3,1} = 'dKroSBL';
supprecovery{3,2} = zeros(JJ,AVG);
supprecovery{4,1} = 'classicSBL';
supprecovery{4,2} = zeros(JJ,AVG);
supprecovery{5,1} = 'OMP';
supprecovery{5,2} = zeros(JJ,AVG);
supprecovery{6,1} = 'dOMP';
supprecovery{6,2} = zeros(JJ,AVG);
supprecovery{7,1} = 'dFISTA';
supprecovery{7,2} = zeros(JJ,AVG);

% symbol error rate
ser = cell(5,2);
ser{1,1} = 'SVD-KroSBL';
ser{1,2} = zeros(JJ,AVG);
ser{2,1} = 'AM-KroSBL';
ser{2,2} = zeros(JJ,AVG);
ser{3,1} = 'dKroSBL';
ser{3,2} = zeros(JJ,AVG);
ser{4,1} = 'classicSBL';
ser{4,2} = zeros(JJ,AVG);
ser{5,1} = 'OMP';
ser{5,2} = zeros(JJ,AVG);
ser{6,1} = 'dOMP';
ser{6,2} = zeros(JJ,AVG);
ser{7,1} = 'dFISTA';
ser{7,2} = zeros(JJ,AVG);

% time
time = cell(5,2);
time{1,1} = 'SVD-KroSBL';
time{1,2} = zeros(JJ,AVG);
time{2,1} = 'AM-KroSBL';
time{2,2} = zeros(JJ,AVG);
time{3,1} = 'dKroSBL';
time{3,2} = zeros(JJ,AVG);
time{4,1} = 'classicSBL';
time{4,2} = zeros(JJ,AVG);
time{5,1} = 'OMP';
time{5,2} = zeros(JJ,AVG);
time{6,1} = 'dOMP';
time{6,2} = zeros(JJ,AVG);
time{7,1} = 'dFISTA';
time{7,2} = zeros(JJ,AVG);

% noise_var
noise_var_error = cell(4,2);
noise_var_error{1,1} = 'SVD-KroSBL';
noise_var_error{1,2} = zeros(JJ,AVG);
noise_var_error{2,1} = 'AM-KroSBL';
noise_var_error{2,2} = zeros(JJ,AVG);
noise_var_error{4,1} = 'classicSBL';
noise_var_error{4,2} = zeros(JJ,AVG);