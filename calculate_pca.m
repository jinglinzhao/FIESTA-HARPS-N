X       = dlmread('X.txt')';
X_err   = dlmread('X_err.txt')';

% W       = 1 ./ X_err;
% W(7,:) = W(7,:)/100;
W = ones(size(W));

ncomp = 11;
niter = 2000;
nrefine = 2000;
xi = 0;

[P, C] = WPCA(X, W, ncomp, niter, nrefine, xi);
% writematrix(C','pca_score1.txt')
% writematrix(C','pca_score2.txt') % W(7,:) = W(7,:)/100;
writematrix(C','pca_score3.txt') % W = ones(size(W));