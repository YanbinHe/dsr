function [metrics] = classicSBL_noiseless(y,A,N,R_max,x)

time = 0;
A_ori = A;
gammac = 1*ones(N^3,1); % initialize the prior
keep_listc = 1:N^3;

normc = 1;
itrc = 1;

r_abso = 1e-5;
thres = 1e-4;

while(itrc < R_max + 1 && normc > thres) % do the iteration
    itrc;
    
    tic;
    
    gammal = zeros(N^3,1);
    gammal(keep_listc) = gammac;
    gamma_old = gammal;


    if min(abs(gammac)) < r_abso %*max(abs(gammac))
        indexc = find(gammac > r_abso); % should keep *max(abs(gammac))
        gammac = gammac(indexc);
        A = A(:,indexc);
        keep_listc = keep_listc(indexc);
    end

    % update the posterior
    Gammac = diag(gammac);
    Gammacsqrt = diag(gammac.^0.5);
    
    temp = Gammacsqrt*pinv(A*Gammacsqrt);
    Sigma_x = Gammac - temp*A*Gammac;
    x_re = temp * y;
    
    % updating gamma
    Sigma_x_diag = abs(real(diag(Sigma_x)));
    gammac = Sigma_x_diag + abs(x_re).^2;
    
    time = time + toc;

    gammal = zeros(N^3,1);
    gammal(keep_listc) = gammac;

    normc = norm(gammal-gamma_old)/norm(gamma_old);
    
    itrc = itrc + 1; 
end

% re-estimate
Gammac = diag(gammal);
Gammacsqrt = diag(gammal.^0.5);

temp = Gammacsqrt*pinv(A_ori*Gammacsqrt);
Sigma_x = Gammac - temp*A_ori*Gammac;
x_rel = temp * y;

errorc = norm(x_rel - x)/norm(x)

srr = recover_rate(x_rel,x);
metrics={'error',errorc;
    'vector',x_rel;
    'support recovery rate',srr;
    'time',time
    };

end

