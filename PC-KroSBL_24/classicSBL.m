function [time,x_rel,noise_var] = classicSBL(noise_var,numItr,Res1,H,y,SNR)

time = 0;

gammac = 1*ones(Res1^3,1); % initialize the prior
keep_listc = 1:Res1^3;
% noise_var = 1e-1;

normc = 1;
itrc = 1;
thres = 1e-4;
if SNR >= 20
    r = 1e-4;
else
    r = 1e-4;
end

while(itrc < numItr + 1 && normc > thres) % do the iteration
    itrc;
    
    tic;
    
    gammacl = zeros(Res1^3,1);
    gammacl(keep_listc) = gammac;
    gammac_old = gammacl;

    if min(abs(gammac)) < r*max(abs(gammac))
         
        indexc = find(gammac > r*max(abs(gammac))); % should keep
        gammac = gammac(indexc);
        H = H(:,indexc);
        keep_listc = keep_listc(indexc);   
    end
    gammad = gammac;

    % update the posterior mean and variance
    Gammac = diag(gammac);
    invPhic = H'* (noise_var * eye(size(H,1)) + H * Gammac * H')^-1 * H;
    Sigma_x = Gammac - Gammac * invPhic * Gammac;
    x_re = noise_var^-1 * Sigma_x * H' * y;

    % updating gamma
    Sigma_x_diag = real(diag(Sigma_x));
    gammac = Sigma_x_diag + abs(x_re).^2;

    % updating noise variance
%     temp = sum(ones(length(gammad),1) - gammad.^-1 .* Sigma_x_diag);
%     noise_var = (norm(y - H*x_re)^2 + noise_var*temp)/length(y);
   
    time = time + toc;

    gammal = zeros(Res1^3,1);
    gammal(keep_listc) = gammac;

    normc = norm(gammal-gammac_old)/norm(gammac_old);
    norm(gammal-gammac_old)/norm(gammac_old);
    itrc = itrc + 1; 
end

% re-estimating
Gammac = diag(gammac);
invPhic = H'* (noise_var * eye(size(H,1)) + H * Gammac * H')^-1 * H;
Sigma_x = Gammac - Gammac * invPhic * Gammac;
x_re = noise_var^-1 * Sigma_x * H' * y;

x_rel = zeros(Res1^3,1);
x_rel(keep_listc) = x_re;
end

