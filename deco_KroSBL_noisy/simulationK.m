function [resultK] = simulationK(A1_ori,A2_ori,A3_ori,M1,m_loc,lenK,N,R_max,K,SNR,epi)

resultK = cell(1,lenK);

k = 1;
supp1 = randsample(1:N, K(k));
supp2 = randsample(1:N, K(k));
supp3 = randsample(1:N, K(k));

x1 = zeros(N,1);
x2 = zeros(N,1);
x3 = zeros(N,1);

for k = 1:lenK
    
    if k == 1
        x1(supp1) = randn(K(k),1);
        x2(supp2) = randn(K(k),1);
        x3(supp3) = randn(K(k),1);
    
    else
        
        zo1 = find(x1 == 0);
        zo2 = find(x2 == 0);
        zo3 = find(x3 == 0);
        
        new1 = randsample(zo1, 1);
        new2 = randsample(zo2, 1);
        new3 = randsample(zo3, 1);
        
        x1(new1) = randn();
        x2(new2) = randn();
        x3(new3) = randn();
    end

    x = kron(x1,kron(x2,x3));
    suppTrue = find(abs(x)>0);

    y1 = A1_ori(1:M1(m_loc),:)*x1;
    y2 = A2_ori*x2;
    y3 = A3_ori*x3;

    y = kron(y1,kron(y2,y3));
    power_symbol = norm(y)^2/length(y);
    noisesq = power_symbol/(10^(SNR/10));
    noise = sqrt(noisesq)*randn(size(y));
    y = y + noise;

    % measurements
    A1 = A1_ori(1:M1(m_loc),:);
    A2 = A2_ori(:,:);
    A3 = A3_ori(:,:);

    A = kron(A1,kron(A2,A3));

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
    result{7} = noisesq;
    resultK{k} = result;
end
end