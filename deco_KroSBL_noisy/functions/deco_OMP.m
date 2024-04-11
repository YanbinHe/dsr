function [metrics] = deco_OMP(y,A1,A2,A3,N,x,epi)

A{1} = A1;
A{2} = A2;
A{3} = A3;

dim1 = size(A1,1);
dim2 = size(A2,1);
dim3 = size(A3,1);

tic;
yten = cell(3,1);
mat1 = reshape(y,dim2*dim3,dim1);
[mat1l,mat1v,mat1r] = svd(mat1.');
yten{1} = mat1l(:,1);
mat2 = reshape(mat1v(1,1)*mat1r(:,1), dim3, dim2);
[mat2l,mat2v,mat2r] = svd(mat2.');
yten{2} = mat2l(:,1);
yten{3} = mat2v(1,1)*mat2r(:,1);
x_re = cell(3,1);
time_deco = toc;

tic;
for i = 1:3
    itr = 0;
    M = size(A{i},1);
    x_re{i} = zeros(N,1);
    r = yten{i};
    Ti = [];
    % norm of each column 
    norm_A = zeros(1,N);
    for j = 1:N
        norm_A(1,j) = norm(A{i}(:,j),2);
    end 
    
    while(norm(r) > epi && itr < M)
        
        % Normalised correlation of each col of A with r
        aj = abs((r'*A{i}) ./ norm_A);
        [ajmax, j] = max(aj);
        Ti = [Ti j];
        At = A{i}(:,Ti);
        s = (At' * At) \ At' * yten{i};
        r = yten{i} - At*s;

        itr = itr + 1;
    end

    x_re{i}(Ti) = s;
end

x_relv = kron(x_re{1},kron(x_re{2},x_re{3}));
time = toc;

errordo = (norm(x_relv - x)/norm(x))^2

srr = recover_rate(x_relv,x);
metrics={'error',errordo;
    'vector',x_relv;
    'support recovery rate',srr;
    'time',(time+time_deco)
    };


end