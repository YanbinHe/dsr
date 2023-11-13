function r = recover_rate(x_re,xt)
% computing thr support recovery support

% using the equation defined in "Alternative to Extended Block Sparse 
% Bayesian Learning and Its Relation to Pattern-Coupled Sparse Bayesian 
% Learning"

delta = 1e-4;

suppxre = find(abs(x_re) > delta); % to eliminate the potential small terms that come from estimation inaccuracy     *max(abs(x_re))
suppxt = find(abs(xt) > 0);

nomi = length(intersect(suppxre,suppxt));
deno = length(union(suppxre,suppxt));% the union of two sets

r = nomi/deno;
end

