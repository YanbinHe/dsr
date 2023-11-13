% noisy compare
lenK = length(K);
lenM = length(M1);

for avg = 1:AVG
    % generate measuring dictionaries
    A1_ori = randn(M1(end),N);
    A2_ori = randn(M2,N);
    A3_ori = randn(M3,N);
    
    %% fig 1/2 (b) nmse and srr vs measurement      
    k = 2;

    x1 = zeros(N,1);
    x2 = zeros(N,1);
    x3 = zeros(N,1);

    supp1 = randsample(1:N, K(k));
    x1(supp1) = randn(K(k),1);
    supp2 = randsample(1:N, K(k));
    x2(supp2) = randn(K(k),1);
    supp3 = randsample(1:N, K(k));
    x3(supp3) = randn(K(k),1);

    x = kron(x1,kron(x2,x3));
    suppTrue = find(abs(x)>0);

    y1 = A1_ori*x1;
    y2 = A2_ori*x2;
    y3 = A3_ori*x3;

    y_ori = kron(y1,kron(y2,y3));

    for m = 1:lenM
        simulationM;
    end

    %% fig 1/2 (c) nmse and srr vs sparsity level    
    m = 2;
   
    for k = 1:lenK
        simulationK;
    end
    
    %% save data for each trial
    filename = ['./results/noise_compare_', num2str(avg),'.mat'];
    save(filename, 'resultS','resultM','resultK')
end
