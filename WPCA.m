function [P, C] = WPCA(X, W, ncomp, niter, nrefine, xi)
    nvar = size(X, 1); nobs = size(X,2); 
    P = zeros(nvar, ncomp); C = zeros(ncomp, nobs); 
    ws = sum(W, 2);
    covar = ((ws * ws').^xi).*((X.*W)*(X.*W)')./(W*W');
    covar(isnan(covar))=0;
    
    for i=1:ncomp
        u = ones(nvar,1) / nvar;
        for j=1:niter
            u = covar * u;
            d = sqrt(u' * u);
            u = u / d;
        end
        d = u'* covar * u; 
        for j = 1:nrefine
%             u = inv(covar - d * eye(nvar)) * u;
            u = (covar - d * eye(nvar)) \ u;
            u = u / sqrt(u' * u);
            d = u' * covar * u; 
        end
        covar = covar - u * d * u';
        P(:,i) = u;
    end

    for i=1:nobs
        w = diag(W(:,i)).^2;
        C(:,i) = (P' * w * P) \ (P' * w * X(:,i));
    end
end

