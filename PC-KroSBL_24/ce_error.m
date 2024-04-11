function [error, H_re] = ce_error(x_rel,Res1,K2,irs_pattern,A_irs_a,A_irs_d,A1,A2,H1,H2)

for i = 1:K2
    Htt = irs_pattern(:,i).'*kr(A_irs_a.',A_irs_d').';
    Ht(:,:,i) = Htt(:,1:Res1);
    H_re(:,:,i) = kron(kron(Ht(:,:,i),conj(A1)),A2)*x_rel;
    Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    nmse(i) = (norm(H_re(:,:,i) - Htrue(:,:,i),'fro')/norm(Htrue(:,:,i),'fro'))^2;
end
error = sum(nmse)/K2
end