function [time,x_relv] = deco_kroSBL(numItr,Res1,A,y)

num_dict = 3;
time = zeros(num_dict,1);
thres = 1e-4;

dim1 = size(A{1},1);
dim2 = size(A{2},1);
dim3 = size(A{3},1);

A_ori = A;
x_re = cell(3,1);
gammad_old = cell(3,1);
Gammac = cell(3,1);
Gammacsqrt = cell(3,1);
Sigma_x = cell(3,1);
Sigma_x_diag = cell(3,1);

gammad{1} = ones(Res1,1); % initialize the prior
gammad{2} = ones(Res1,1);
gammad{3} = ones(Res1,1);

keep_listc{1} = 1:Res1;
keep_listc{2} = 1:Res1;
keep_listc{3} = 1:Res1;

noise_var = cell(3,1);
noise_var{1} = 10^-(0.33);
noise_var{2} = 10^-(0.33);
noise_var{3} = 10^-(0.33);
% decompose noisy y
% as the model of noise after decomposition is unknown, we treat this as a
% problem with noise estimating
tic;
yten = cell(3,1);
mat1 = reshape(y,dim2*dim3,dim1);
[mat1l,mat1v,mat1r] = svd(mat1.');
yten{1} = mat1l(:,1);
mat2 = reshape(mat1v(1,1)*conj(mat1r(:,1)), dim3, dim2);
[mat2l,mat2v,mat2r] = svd(mat2.');
yten{2} = mat2l(:,1);
yten{3} = mat2v(1,1)*conj(mat2r(:,1));
time_deco = toc;

% control the convergence status
normc = ones(num_dict,1);
itrc = ones(num_dict,1);
r_rele{1} = 1e-3;
r_rele{2} = 1e-3;
r_rele{3} = 1e-3;

for sys_idx = 1:num_dict
    while(itrc(sys_idx) < numItr + 1 && normc(sys_idx) > thres) % do the iteration

        tic;
        
        
        gammal = zeros(Res1^3,1);
        gammal(keep_listc{sys_idx}) = gammad{sys_idx};
        gammad_old{sys_idx} = gammal;

        if min(abs(gammad{sys_idx})) < r_rele{sys_idx}*max(abs(gammad{sys_idx}))
            indexc = find(gammad{sys_idx} > r_rele{sys_idx}*max(abs(gammad{sys_idx}))); % should keep 
            gammad{sys_idx} = gammad{sys_idx}(indexc);
            A{sys_idx} = A{sys_idx}(:,indexc);
            keep_listc{sys_idx} = keep_listc{sys_idx}(indexc);
        end
        gammaold = gammad;
        % update the posterior
        Gammac{sys_idx} = diag(gammad{sys_idx});
        invPhic = A{sys_idx}'* (noise_var{sys_idx} * eye(size(A{sys_idx},1)) + A{sys_idx} * Gammac{sys_idx} * A{sys_idx}')^-1 * A{sys_idx};
        Sigma_x{sys_idx} = Gammac{sys_idx} - Gammac{sys_idx} * invPhic * Gammac{sys_idx};
        x_re{sys_idx} = noise_var{sys_idx}^-1 * Sigma_x{sys_idx} * A{sys_idx}' * yten{sys_idx};
        
        % updating gamma 
        Sigma_x_diag{sys_idx} = real(diag(Sigma_x{sys_idx}));
        gammad{sys_idx} = Sigma_x_diag{sys_idx} + abs(x_re{sys_idx}).^2;
        
        % updating noise variance
        temp = sum(ones(length(gammaold{sys_idx}),1) - gammaold{sys_idx}.^-1 .* Sigma_x_diag{sys_idx});
        noise_var{sys_idx} = (norm(yten{sys_idx} - A{sys_idx}*x_re{sys_idx})^2 + noise_var{sys_idx}*temp)/length(yten{sys_idx});
        

        time(sys_idx) = time(sys_idx) + toc;

        gammal = zeros(Res1^3,1);
        gammal(keep_listc{sys_idx}) = gammad{sys_idx};

        normc(sys_idx) = norm(gammal-gammad_old{sys_idx})/norm(gammad_old{sys_idx});
        norm(gammal-gammad_old{sys_idx})/norm(gammad_old{sys_idx});
        itrc(sys_idx) = itrc(sys_idx) + 1; 
    end
end

x_rel = cell(3,1);
% re-estimate
for sys_idx = 1:num_dict
    Gammac{sys_idx} = diag(gammad{sys_idx});
    invPhic = A{sys_idx}'* (noise_var{sys_idx} * eye(size(A{sys_idx},1)) + A{sys_idx} * Gammac{sys_idx} * A{sys_idx}')^-1 * A{sys_idx};
    Sigma_x{sys_idx} = Gammac{sys_idx} - Gammac{sys_idx} * invPhic * Gammac{sys_idx};
    x_re{sys_idx} = noise_var{sys_idx}^-1 * Sigma_x{sys_idx} * A{sys_idx}' * yten{sys_idx};
    
    x_rel{sys_idx} = zeros(Res1,1);
    x_rel{sys_idx}(keep_listc{sys_idx}) = x_re{sys_idx};
end

x_relv = kron(x_rel{1},kron(x_rel{2},x_rel{3}));
time = sum(time) + time_deco;

end


