function [metrics] = deco_fista_lasso(y,A1,A2,A3,N,x)

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
time_deco = toc;

results = cell(3,1);

for i = 1:3
    results{i} = fista_lasso(A{i},yten{i},1);
end

x_relv = kron(results{1}{1,2},kron(results{2}{1,2},results{3}{1,2}));

errordf = norm(x_relv - x)/norm(x)

srr = recover_rate(x_relv,x);
metrics={'error',errordf;
    'vector',x_relv;
    'support recovery rate',srr;
    'time',(results{1}{2,2}+results{2}{2,2}+results{3}{2,2}+time_deco)
    };



end