function [time,x_rel,noise_var] = kroSBL_am(noise_var,Res1,numItr,H,H_p1,H_p2,A2,y,SNR)

time = 0;

A2_ori = A2; % a copy of the dictionary
thres_inner = 1e-5;
thres = 1e-4;

gamma1 = 1*ones(Res1,1); % prior on the gamma1
N1 = norm(gamma1);
gamma1 = gamma1/N1;
gamma2 = 1*ones(Res1,1); % prior on the gamma2
N2 = norm(gamma2);
gamma2 = gamma2/N2;
gamma3 = N1*N2*1*ones(Res1,1); % prior on the gamma3
gamma = kron(kron(gamma1,gamma2),gamma3);

keep_list = 1:Res1^3;
keep_list1 = 1:Res1;
keep_list2 = 1:Res1;% absolute position
keep_list3 = 1:Res1;% absolute position
% noise_var = 1e-1;

norm1 = 1;
itr2 = 1;

r_absolute = 1e-4;
if SNR >= 25
    r1 = 1e-3;
else
    r1 = 1e-2;
end
r2 = 1e-1;
r3 = 1e-3;

while(itr2 < numItr + 1 && norm1 > thres) % do the iteration
    tic;
    gammare = zeros(Res1^3,1);
    gammare(keep_list) = gamma;
    gamma_old = gammare;

    if min(abs(gamma1)) < r_absolute || min(abs(gamma2)) < r_absolute || min(abs(gamma3)) < r_absolute || min(abs(gamma1)) < r1*max(abs(gamma1)) || min(abs(gamma2)) < r2*max(abs(gamma2)) || min(abs(gamma3)) < r3*max(abs(gamma3)) 
        % if the first condition happens, it means that some block has to
        % be zero, if the second happens, it means that some entries in
        % each block should be zero.
        
        l1 = length(gamma1);
        l2 = length(gamma2);
        l3 = length(gamma3);   
        
        index1_dele = find(gamma1 < r1*max(abs(gamma1)));
        index1_sub = find(gamma1 < r_absolute);
        index1_dele = sort(unique([index1_dele;index1_sub]));
        
        index1 = 1:l1;
        index1(index1_dele) = [];% should keep relative position
        
        
        index2_dele = find(gamma2 < r2*max(abs(gamma2))); 
        index2_sub = find(gamma2 < r_absolute);
        index2_dele = sort(unique([index2_dele;index2_sub]));
        
        index2 = 1:l2;
        index2(index2_dele) = [];% should keep relative position
        
        index3_dele = find(gamma3 < r3*max(abs(gamma3))); 
        index3_sub = find(gamma3 < r_absolute);
        index3_dele = sort(unique([index3_dele;index3_sub]));
        
        index3 = 1:l3;
        index3(index3_dele) = [];% should keep relative position
        
        index = [];%zeros(l1*l2,1);
        
        for i =1:length(index1)
            for j = 1:length(index2)
                for k = 1:length(index3)
                    index = [index;(((index1(i)-1)*l2+index2(j)-1)*l3+index3(k))];
                end
            end
        end
        
        H = H(:,index);
        H_p1 = H_p1(:,index1);
        H_p2 = H_p2(:,index2);
        A2 = A2(:,index3);
        
        keep_list = keep_list(index);
        keep_list1 = keep_list1(index1);
        keep_list2 = keep_list2(index2);
        keep_list3 = keep_list3(index3);

        % prune gamma and corresponding entries in Sigma and mu
        gamma1 = gamma1(index1);
        gamma2 = gamma2(index2);
        gamma3 = gamma3(index3);

    end
    gammad = kron(gamma1,kron(gamma2,gamma3));
    
    % update the posterior mean and variance
    % compute the posterior
    [x_re,Sigma_x] = posterior_compute2(noise_var,gamma1,gamma2,gamma3,H_p1,H_p2,A2,H,y);
    
    Lambda = real(diag(Sigma_x)) + abs(x_re).^2;  

    normi1 = 1;
    l1 = size(gamma1,1);
    l2 = size(gamma2,1);
    l3 = size(gamma3,1);
    
    while(normi1 > thres_inner)
        
        gamma_old_inner = kron(gamma1,kron(gamma2,gamma3));
        
        % update
        gamma1 = (l2*l3)^-1*(kron(eye(l1),kron(gamma2.^-1,gamma3.^-1)))'*Lambda;

        gamma2 = (l1*l3)^-1*(kron(gamma1.^-1,kron(eye(l2),gamma3.^-1)))'*Lambda;

        gamma3 = (l1*l2)^-1*(kron(gamma1.^-1,kron(gamma2.^-1,eye(l3))))'*Lambda;
        
        % projection
        n1 = norm(gamma1);           
        gamma1 = gamma1/n1;     
        n2 = norm(gamma2);
        gamma2 = gamma2/n2;
        gamma3 = gamma3*n1*n2;
        
        normi1 = norm(kron(kron(gamma1,gamma2),gamma3) - gamma_old_inner)/norm(gamma_old_inner);
    end

    % update noise variance
%     temp = sum(ones(length(gammad),1) - gammad.^-1 .* real(diag(Sigma_x)));
%     noise_var = (norm(y - H*x_re)^2 + noise_var*temp)/length(y);
    
    time = time + toc;   
    
    gamma = kron(gamma1,kron(gamma2,gamma3));
    gammare = zeros(Res1^3,1);
    gammare(keep_list) = gamma;  
    
    norm1 = norm(gammare - gamma_old)/norm(gamma_old);
    norm(gammare - gamma_old)/norm(gamma_old);
    
    itr2 = itr2 + 1;
end
% re-estimate
[x_re,~] = posterior_compute2(noise_var,gamma1,gamma2,gamma3,H_p1,H_p2,A2,H,y);

x_rel = zeros(Res1^3,1);
x_rel(keep_list) = x_re;
end

