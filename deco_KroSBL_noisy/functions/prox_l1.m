function N = prox_l1(u, alpha)
    N = length(u);
    re = zeros(N,1);
    upos = find(abs(u)>=alpha); 
    output = exp(1i * angle(u)).*(abs(u) - alpha);
    re(upos) = output(upos);
end