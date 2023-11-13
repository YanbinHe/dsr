function [time,x_re] = OMP(Res1,H,y,klevel)

time = 0;

tic;

residual = y;
set = [];
x_re = zeros(Res1^3,1);

while (length(set) < klevel + 1 && norm(residual)^2/length(y) > 1e-3)

    inner_pro = abs(H'*residual);
    [~,in] = max(inner_pro);
    set = [set in];
    
    Hset = H(:,set);
    xset = pinv(Hset)*y;

    x_re(set) = xset;
    residual = y - H*x_re;
end

time = time + toc;
end

