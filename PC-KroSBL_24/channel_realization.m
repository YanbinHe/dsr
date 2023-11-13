% generate channel realization
AoD_ms = randsample([1:Res1], 1);
g1 = zeros(Res1,1);
g1(AoD_ms) = 1;

Smpl = randsample([0:Res1-angle_spread], 1);
AoA_irs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;
ga = zeros(Res1,1);
ga(AoA_irs) = 1;

AoD_irs = randsample([1:Res1], 1);
gd = zeros(Res1,1);
gd(AoD_irs) = 1;

Smpl = randsample([0:Res1-angle_spread], 1);
AoA_bs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;  
g2 = zeros(Res1,1);
g2(AoA_bs) = 1;

% generate the true support, for the support recovery computation
gi = zeros(Res1,1);
indexi = [(AoA_irs(1)-AoD_irs(end)):1:(AoA_irs(end)-AoD_irs(1))] + Res1;
indexmodi = find(indexi > Res1);
indexi(indexmodi) = indexi(indexmodi) - Res1;
gi(indexi) = 1; 
gi = (fliplr(gi'))';

% suppTrue is a vector where the true position has value 1
suppTrue = kron(gi,kron(g1,g2));
%%
% dictionary generation
A2 = generate_dict(bs_ante,Res1);
A_irs_d = generate_dict(irs_ele,Res1);
A_irs_a = generate_dict(irs_ele,Res1);
A1 = generate_dict(ms_ante,Res1); 

% path gain CN(0,1)
alpha_aco = (1*randn(angle_spread,1) + 1i*randn(angle_spread,1))/sqrt(2); % from ms to irs
alpha_a = zeros(Res1,1);
alpha_a(AoA_irs) = alpha_aco;

alpha_2co = (1*randn(angle_spread,1) + 1i*randn(angle_spread,1))/sqrt(2); % from irs to bs
alpha_2 = zeros(Res1,1);
alpha_2(AoA_bs) = alpha_2co;

H2 = sqrt(bs_ante*irs_ele/angle_spread)*A2*(g2.*alpha_2)*gd'*A_irs_d';
H1 = sqrt(ms_ante*irs_ele/angle_spread)*A_irs_a*(ga.*alpha_a)*g1'*A1';

% make sure that somehow the irs should point to the receiver, otherwise
% the received signal is too weak. But this is only for the simulation.
% This is not the case when it comes to real world.1/sqrt(irs_ele)*
for i = 1:max(K2)
    Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    while norm(Htrue(:,:,i),'fro')<1e-5
        irs_pattern(:,i) =  1/sqrt(irs_ele)*((rand(irs_ele,1) > 0.5)*2-1); 
        Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    end
end

% generate signal
y_bar = kr(X.'*H1.',H2)*irs_pattern;
signal_power = norm(vec(y_bar))^2/length(vec(y_bar)); %average signal power per asymbol
