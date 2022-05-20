function X = computeX(Xinit, V, b, N, I, NT, sigma, dt)

X = zeros(N*I, NT); % Note: we store the states resulting from different 
                    % initial conditions (i.e. input data)
                    % in one long vector of length N*I. 
X(:,1) = Xinit;
for kk = 1:(NT-1)
    for ii = 1:I
        ind = N*(ii-1)+(1:N);
        X(ind,kk+1) = X(ind,kk) + TODO;
    end
end