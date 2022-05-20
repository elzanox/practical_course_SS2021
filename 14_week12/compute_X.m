function [X, duration] = compute_X(E, A, Xinit, B, U, NT, dt)

X = zeros(length(Xinit), NT);
X(:,1) = Xinit;
for ii = 2:NT
    X(:,ii) = TODO;
end