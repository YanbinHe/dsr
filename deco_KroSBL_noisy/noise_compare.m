% noise compare
lenK = length(K);
lenM = length(M1);
lenS = length(SNR);

for avg = 1:AVG
    resultM = cell(1,lenM);
    resultS = cell(1,lenS);
    
    % generate measuring dictionaries
    A1_ori = randn(M1(end),N);
    A2_ori = randn(M2,N);
    A3_ori = randn(M3,N);

    %% nmse and srr vs snr
    k_loc = 2;
    m_loc = 6; % 12 measurement level

    x1 = zeros(N,1);
    x2 = zeros(N,1);
    x3 = zeros(N,1);
    
    % generate sparse vectors
    supp1 = randsample(1:N, K(k_loc));
    x1(supp1) = randn(K(k_loc),1);
    supp2 = randsample(1:N, K(k_loc));
    x2(supp2) = randn(K(k_loc),1);
    supp3 = randsample(1:N, K(k_loc));
    x3(supp3) = randn(K(k_loc),1);

    A1 = A1_ori(1:M1(m_loc),:);
    A2 = A2_ori;
    A3 = A3_ori;
    
    y1 = A1*x1;
    y2 = A2*x2;
    y3 = A3*x3;

    x = kron(x1,kron(x2,x3));
    y = kron(y1,kron(y2,y3));
    A = kron(A1,kron(A2,A3));
    suppTrue = find(abs(x)>0);
    %%
    for s = 1:lenS
        resultS{s} = simulationS(y,A1,A2,A3,A,N,R_max,x,SNR(s),epi);
    end
    
    %% nmse and srr vs measurement 
    SNR_loc = 20;
    k_loc = 2; % 3 non-zeros in each vector

    x1 = zeros(N,1);
    x2 = zeros(N,1);
    x3 = zeros(N,1);
    
    % generate sparse vectors
    supp1 = randsample(1:N, K(k_loc));
    x1(supp1) = randn(K(k_loc),1);
    supp2 = randsample(1:N, K(k_loc));
    x2(supp2) = randn(K(k_loc),1);
    supp3 = randsample(1:N, K(k_loc));
    x3(supp3) = randn(K(k_loc),1);
    
    y1 = A1_ori*x1;
    y2 = A2_ori*x2;
    y3 = A3_ori*x3;

    x = kron(x1,kron(x2,x3));
    suppTrue = find(abs(x)>0);
    for m = 1:lenM
        resultM{m} = simulationM(y1,y2,y3,A1_ori,A2_ori,A3_ori,M1,m,N,R_max,x,SNR_loc,epi);
    end

    %% fig 1/2 (c) nmse and srr vs sparsity level    
    SNR_loc = 20;
    m_loc = 6; % 12 measurement level
   
    resultK = simulationK(A1_ori,A2_ori,A3_ori,M1,m_loc,lenK,N,R_max,K,SNR_loc,epi);

    %% save data for each trial
    fname = ['./results/noisy_compare_', num2str(avg),'.mat'];
    save(fname,"resultK","resultM","resultS")
end
