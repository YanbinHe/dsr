function [time,x_relv] = deco_OMP(Res1,A,y,epi)

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
x_re = cell(3,1);
time_deco = toc;

time = zeros(3,1);

for i = 1:3
    [time(i), x_re{i}] = omp2(A{i},yten{i},epi,size(A{i},1));
end

x_relv = kron(x_re{1},kron(x_re{2},x_re{3}));
time = sum(time) + time_deco;
end