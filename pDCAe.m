function [gg,out] = pDCAe(Y, X, lambda, K, pars)
% A solver for TgSLOPE problem:
%    min 1/2||Y_3-XGH'||_F^2 + sum_j (lambda_j * ||G||_[j]),  s.t. H'H=I,
%    Y: n*p1*p2, X: n*p, G: p*K, H: p1p2*K
% Inputs:
%     Y:        Responses
%     X:        Design matrix
%     lambda:   Regularization parameter vector of TgSLOPE
%     K:        CP rank,<= min(p,p1*p2)
%     pars:     Parameters are all OPTIONAL
%               pars.iteron --  =1. Results will  be shown for each iteration (default)
%                               =0. Results won't be shown for each iteration
%               pars.maxit  --  Maximum nonumber of iteration.  (default 1000)
%               pars.tol    --  Tolerance of stopping criteria. (default 1e-5)
%
% Outputs:
%     gg:           The row L2-norm vector of solution G,(||G(1,:)||,...||G(p,:)||)'
%     out.G:        The solution G
%     out.H:        The solution H
%     out.runtime   CPU time

% -------------------------------------------------------------
% Start timer
% -------------------------------------------------------------
t0 = tic();

% -------------------------------------------------------------
% Parse parameters
% -------------------------------------------------------------
warning off;
if nargin<4; error('Imputs are not enough!\n'); end
if nargin<5; pars=struct(); end
if isfield(pars,'iteron');iteron = pars.iteron; else; iteron = 1;        end
if isfield(pars,'maxit'); maxit  = pars.maxit;  else; maxit  = 5e4;      end
if isfield(pars,'tol');   tol    = pars.tol;    else; tol = 1e-6;        end
if isfield(pars,'ginit'); GInit  = pars.GInit;  else; GInit = [];        end
if isfield(pars,'hinit'); Q1Init  = pars.HInit;  else; Q1Init = [];      end

[~,p]     = size(X);
X2        = X'*X;
YX        = Y'*X;
if p<= 1000
    L         = norm(X2,2);
else
    opts.issym = 1;
    L = eigs(X2,1,'LM',opts);
end
if (isempty(GInit)), GInit = ones(p,K); end
if (isempty(Q1Init))
    [U,S,V]= svd(YX*GInit,'econ');
    Q1Init = YX'*U*V';
end

% Initialization
G       = GInit;
Q1      = Q1Init;
D       = GInit;
f       = objval(lambda,X,G,S);
beta    = 1;   % beta = 0: without extrapolation
if iteron
    fprintf('\n Iter      Obj       Residual       Time\n');
    fprintf('----------------------------------------  \n');
end

% Main loop
for iter = 1: maxit
    G_old = G;
    Q1_old = Q1;
    beta_old = beta;
    D_old  = D;
    Q_old = L*D_old-X2*D_old+Q1_old;
    QL_old = norm_20(Q_old)/L;
    
    % Proximal operator calculation of the group SLOPE penalty
    g = proxSortedL1(QL_old,lambda/L);
    for j=1:p
        if g(j) ==0
            G(j,:) =  0;
        else
            G(j,:) = g(j)*Q_old(j,:)/norm(Q_old(j,:),2);
        end
    end
    
    [U,S,V] = svd(YX*G,'econ');
    Q1  = YX'*U*V';
    beta = (1+sqrt(1+4*beta_old^2))/2;
    D    = G+((beta_old-1)/beta)*(G-G_old);
    f    = objval(lambda,X,G,S);
    
    % Stop criteria
    residual = norm(G-G_old,'fro')/max(norm(G_old,'fro'),1);
    if iteron && mod(iter,50)==0
        fprintf('%4d       %5.2e     %5.2e      %5.2fsec\n',iter,f,residual,toc(t0));
    end
    if residual<tol
        fprintf('%4d      %5.2e    %5.2e     %5.2fsec\n',iter,f,residual,toc(t0));
        break;
    end
end
H  = U*V';
gg = g;

% Information structure
if iteron
    fprintf('--------------------------------------------------\n');
    out.G         = G;
    out.H         = H;
    out.B         = G*H';
    out.runtime   = toc(t0);
end
end

function obj = objval(lambda,X,G,S)

% Compute objective value

obj = norm(X*G,'fro')^2/2+lambda*norm_20(G)-trace(S);
end

function x = proxSortedL1(y,lambda)

% Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes

% This file is part of SLOPE Toolbox version 1.0.
%
%    The SLOPE Toolbox is free software: you can redistribute it
%    and/or  modify it under the terms of the GNU General Public License
%    as published by the Free Software Foundation, either version 3 of
%    the License, or (at your option) any later version.
%
%    The SLOPE Toolbox is distributed in the hope that it will
%    be useful, but WITHOUT ANY WARRANTY; without even the implied
%    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%    See the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with the SLOPE Toolbox. If not, see
%    <http://www.gnu.org/licenses/>.

% Normalization
lambda = lambda(:);
y       = y(:);
sgn     = sign(y); % Returns phase for complex numbers
[y,idx] = sort(abs(y),'descend');

% Simplify the problem
k = find(y > lambda,1,'last');

% Compute solution and re-normalize
n = numel(y);
x = zeros(n,1);

if (~isempty(k))
    v1 = y(1:k);
    v2 = lambda(1:k);
    v = proxSortedL1Mex(v1,v2);
    x(idx(1:k)) = v;
end

% Restore signs
x = sgn .* x;
end

function b=norm_20(B)

% Compute the L2-norm of the rows of a matrix and reture a vector
p=size(B,1);
z=zeros(p,1);
for i=1:p
    bi=norm(B(i,:),2);
    z(i)=bi;
end
b=z;
end