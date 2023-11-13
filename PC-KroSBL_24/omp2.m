function [time, theta] = omp2(A,y,epi,klevel)
    N = size(A,2);
    M = size(A,1);
    theta = zeros(N,1);
    y = y(:);
    r = y;
    Ti = [];
    
    itr = 0;

    tic;

    % norm of each column 
    norm_A = zeros(1,N);
    for j = 1:N
        norm_A(1,j) = norm(A(:,j),2);
    end 
    
    while(norm(r) > epi && itr < klevel)
        
        % Normalised correlation of each col of A with r
        aj = abs((r'*A) ./ norm_A);
        [ajmax, j] = max(aj);
        Ti = [Ti j];
        At = A(:,Ti);
        s = (At' * At) \ At' * y;
        r = y - At*s;

        itr = itr + 1;
    end
    
    theta(Ti) = s;

    time = toc;
end