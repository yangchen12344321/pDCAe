% Demon of pDCAe to solving TgSLOPE problem with randomly generated data
clc; clear; close all;
n  = 2000;
p  = 1000;
p1 = 10;
p2 = 10;
q  = p1*p2;
K  = 20;
S = [25,50,75,100,125,150,175,200,225,250];
%%
iter     = 100;
alpha1   = 0.05;
alpha2   = 0.1;
nominalFDR_005  = alpha1*(p-S)/p;
nominalFDR_01   = alpha2*(p-S)/p;
Lgths         = K*ones(1,p);
w             = ones(1,p);
FDR_005       = zeros(1,length(S));
FDR_01        = zeros(1,length(S));
stdFDR_005    = zeros(1,length(S));
stdFDR_01     = zeros(1,length(S));
POWER_005     = zeros(1,length(S));
POWER_01      = zeros(1,length(S));
stdPOWER_005  = zeros(1,length(S));
stdPOWER_01   = zeros(1,length(S));
% Signal Strength
aa       = sqrt(4*log(p)/(1-p^(-2/K))-K);
pars.tol = 1e-5;
% lambda
lambda_005  =Lambda(n, alpha1, Lgths, w, 'orthogonal'); % choose type 'orthogonal'or 'gaussian'
lambda_01   = Lambda(n, alpha2, Lgths, w, 'orthogonal');
% main loop for all sparsities
for k = 1:length(S)
    fprintf('sparsity:=%d\n',S(k));
    s          = S(k);
    fdr_005    = zeros(1,iter);
    fdr_01     = zeros(1,iter);
    power_005  = zeros(1,iter);
    power_01   = zeros(1,iter);
    for i = 1:iter
        fprintf('iter:=%d\n',i);
        G         = zeros(p,K);
        I         = randperm(p);
        TG        = I(1:s);
        for j = 1:s
            signals     = abs(rand(1,K)+0.1);
            G(TG(j),:)  = signals*aa*sqrt(K)/norm(signals,2);
        end
        H         = randn(q,q);
        [H,~,~]   = svd(H);
        H         = H(:,1:K);
        %%% Orthogonal design
        X         = randn(n,p);
        X = orth(X);
        %%% Gaussian random design
        %         SIGMA     = zeros(p,p);
        %         for i=1:p
        %             for j=1:p
        %                 SIGMA(i,j) = 0.5^abs(i-j);
        %             end
        %         end
        %         MU      = zeros(1,p);
        %         X       = mvnrnd(MU,SIGMA,n);
        
        Y         = X*G*H'+ randn(n,q);
        gg_005    = pDCAe(Y, X, lambda_005, K, pars);
        gg_005(abs(gg_005)<1e-8) = 0;
        
        gg_01     = pDCAe(Y, X, lambda_01, K, pars);
        gg_01(abs(gg_01)<1e-8)   = 0;
        
        discoveries_005      = find(abs(gg_005)>0);
        discoveries_01       = find(abs(gg_01)>0);
        numb_disc_005        = length(discoveries_005);
        numb_disc_01         = length(discoveries_01);
        true_disc_005        = intersect(discoveries_005,TG);
        true_disc_01         = intersect(discoveries_01,TG);
        false_disc_005       = setdiff(discoveries_005,TG);
        false_disc_01        = setdiff(discoveries_01,TG);
        numb_false_disc_005  = length(false_disc_005);
        numb_false_disc_01   = length(false_disc_01);
        
        fdr_005(i)    = numb_false_disc_005/max(numb_disc_005,1);
        fdr_01(i)     = numb_false_disc_01/max(numb_disc_01,1);
        power_005(i)  = (numb_disc_005-numb_false_disc_005)/s;
        power_01(i)   = (numb_disc_01-numb_false_disc_01)/s;
    end
    FDR_005(k)      = mean(fdr_005);
    FDR_01(k)       = mean(fdr_01);
    stdFDR_005(k)   = std(fdr_005);
    stdFDR_01(k)    = std(fdr_01);
    POWER_005(k)    = mean(power_005);
    POWER_01(k)     = mean(power_01);
    stdPOWER_005(k) = std(power_005);
    stdPOWER_01(k)  = std(power_01);
end
%% saving
cd('...\TgSLOPE\Fig4\Fig4_data\Orth') % Position of File 'Fig4_data'
save lambda_005.txt -ascii lambda_005;
save lambda_01.txt -ascii lambda_01;
save FDR_005.txt -ascii FDR_005;
save FDR_01.txt -ascii FDR_01;
save stdFDR_005.txt -ascii stdFDR_005;
save stdFDR_01.txt -ascii stdFDR_01;
save nominalFDR_005.txt -ascii nominalFDR_005;
save nominalFDR_01.txt -ascii nominalFDR_01;
save POWER_005.txt -ascii POWER_005;
save POWER_01.txt -ascii POWER_01;
save stdPOWER_005.txt -ascii stdPOWER_005;
save stdPOWER_01.txt -ascii stdPOWER_01;