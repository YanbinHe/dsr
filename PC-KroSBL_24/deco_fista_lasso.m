function [time,x_relv] = deco_fista_lasso(A,y)

dim1 = size(A{1},1);
dim2 = size(A{2},1);
dim3 = size(A{3},1);

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

results = cell(3,1);

for i = 1:3
    results{i} = fista_lasso(A{i},yten{i},0.01);
end

x_relv = kron(results{1}{1,2},kron(results{2}{1,2},results{3}{1,2}));
time = results{1}{2,2}+results{2}{2,2}+results{3}{2,2}+time_deco;
end