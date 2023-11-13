function r = recover_rate(x,xtilde)
% computing thr support recovery support
% x: true vector
% xtilde: the estimated one

delta = 1e-4;

suppx= find(abs(x) > 0);
suppxt= find(abs(xtilde) > delta*max(abs(xtilde)));

nomi = length(intersect(suppx,suppxt));
deno = length(union(suppx,suppxt));% the union of two sets

r = nomi/deno;
end

