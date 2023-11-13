function [mean,var_diag] = posterior_compute3(sigmasq,gamma1,gamma2,gamma3,A,B,C,H,y)
% 

gA = (gamma1.*A');
gB = (gamma2.*B');
gC = (gamma3.*C');

G = kron(gamma1,kron(gamma2,gamma3));

Abar = A*gA;
Bbar = B*gB;
Cbar = C*gC;

[U1,P1] = eig(Abar);
[U2,P2] = eig(Bbar);
[U3,P3] = eig(Cbar);

ghhu = kron((gA*U1),kron(gB*U2,gC*U3));
P = kron(diag(P1),kron(diag(P2),diag(P3)));
diaginv = (sigmasq*ones(size(P,1),1)+P).^-1;
ginv = (ghhu*(diaginv.*ghhu'));

hy = H'*y;
var_diag = real(G - diag(ginv));
mean = sigmasq^-1*(G.*hy - (ginv*hy));
end

