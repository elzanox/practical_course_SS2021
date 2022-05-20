function Phi = computePhi(Yout, X, V, b, N, I, NT, dsigma, weights, dt)

Phi = zeros(N*I, NT-1);
Phi(:,end) = TODO;
for kk = NT-2:-1:1
    for ii = 1:I
        ind = N*(ii-1)+(1:N);
        Phi(ind,kk) = Phi(ind,kk+1) + TODO;
    end
end