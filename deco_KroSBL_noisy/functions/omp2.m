function [metrics] = omp2(y,A,x,epi)
    N = size(A,2);
    M = size(A,1);
    x_re = zeros(N,1);
    y = y(:);
    r = y;
    Ti = [];
    
    itr = 0;
    time = 0;
    nnz = length(find(x~=0));

    tic;

    % norm of each column 
    norm_A = zeros(1,N);
    for j = 1:N
        norm_A(1,j) = norm(A(:,j),2);
    end 
    
    while(norm(r) > epi && itr < 4*nnz)
        
        % Normalised correlation of each col of A with r 
        aj = abs((r'*A) ./ norm_A);
        [ajmax, j] = max(aj);
        Ti = [Ti j];
        At = A(:,Ti);
        s = (At' * At) \ At' * y;
        r = y - At*s;

        itr = itr + 1;
    end
    
    x_re(Ti) = s;

    time = time + toc;
    erroro = (norm(x_re - x)/norm(x))^2
    
    srr = recover_rate(x_re,x);
    metrics={'error',erroro;
             'vector',x_re;
             'support recovery rate',srr;
             'time',time
             };
end