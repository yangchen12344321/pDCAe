function [lambda, kink] = Lambda(n, q, Lgths, w, lam_type)

% Lambda sequence for TgSLOPE. 
% n: number of rows of X
% q: target FDR level
% Lghts: vector of groups sizes
% w: vector of groups weights

if nargin ==4
    lam_type = 'gaussian';
end

% Objects
Lgths            = ( Lgths(:) )';
m                = length(Lgths);
w                = ( w(:) )';
lambda           = zeros(1,m);
critical_pvalues = (1:m)*q/m;

% Finding the first lambda
y1        = 1-critical_pvalues(1);
SCALE1    = w.^-1;
lambda(1) = max( SCALE1.*sqrt(chi2inv(y1, Lgths)) );

% Finding rest lambdas
endd = 0;
ii   = 2;
while endd == 0
    y          = 1-critical_pvalues(ii);
    S          = ii-1;
    lam_cur    = lambda(1:(ii-1));
    if strcmp(lam_type, 'orthogonal')
        SCALE = w.^-1;
    else
        SCALE   = w.^-1.*sqrt( (n-Lgths*S)/n+w.^2*(norm(lam_cur,2).^2)/(n-Lgths*S-1) );
    end
    lambda(ii) = max( SCALE.*sqrt(chi2inv(y, Lgths)) );
    if lambda(ii) > lambda(ii-1) || ii>=m
        endd = 1;
    end
    ii = ii + 1;
end

lambda((ii-1):end) = lambda(ii-2);
kink = ii - 2;

end
