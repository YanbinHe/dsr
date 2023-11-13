function metrics = fista_lasso(A,b,lam)
%% initialize vars
time = 0;
AtA = A'*A; Atb = A'*b;
L = norm(AtA, 2)/lam;

x = zeros(size(A, 2), 1);
y = zeros(size(A, 2), 1);

delta = 10;
epsilon = 10^(-2);
num_itrs = 1;

t(num_itrs) = 1;
tic;
%% FISTA
while(delta > epsilon)

    val(num_itrs) = eval_lasso(A, b, lam, x(:, num_itrs));

    u = y - (AtA*x(:, num_itrs) - Atb)/(lam*L);

    x(:, num_itrs + 1) = prox_l1(u, 1/L);

    t(num_itrs+1) = (1 + sqrt(1 + 4*t(num_itrs)^2))/2;

    y = x(:,num_itrs +1) + (1/t(num_itrs+1))*(t(num_itrs) - 1)*(x(:,num_itrs+1) - x(:,num_itrs));

    delta_vec(num_itrs+1) = norm( x(:, num_itrs+1) - x(:, num_itrs), 2);

    delta = delta_vec(num_itrs+1);

    num_itrs = num_itrs + 1;

end

time = time + toc;

x_re = x(:,end);

metrics={
    'vector',x_re;
    'time',time
    };
end