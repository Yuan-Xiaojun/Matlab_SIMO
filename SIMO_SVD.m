function results = SIMO_SVD( optIn, PrioriIn )

%%
M = optIn.M;
N = optIn.N;
L = optIn.L;


A = PrioriIn.A;
s = PrioriIn.s;
rho = PrioriIn.rho ;
X = PrioriIn.X;
Hd = PrioriIn.Hd;
Y = PrioriIn.Y;
Sam = PrioriIn.Sam;
State = PrioriIn.State;
nuw = PrioriIn.noiseVar ;
%%
shat.nm = ones(M,N);
for m = 1 : M
    p = randperm(N) ;
    shat.nm(m , p(1: floor(N * (1 - rho)))) = 0;
end
svar.nm = (rho - rho^2)*ones(M, N);

%% SVD estimate xhat
% [~, ~, Vy] = svd(Y) ;
% xhat = Vy(:,1)' ;
% phase = Sam(1)/xhat(1) ;
% xhat = phase * xhat ;
% X_sam = zeros(length(Sam), L) ;
% for i = 1 : length(Sam)
%        X_sam(i,:) = abs( Sam(i) - xhat );
% end
% [~,I] = min (X_sam);
% xhat = Sam(I);
%  [~ , Xerr] = symerr(X , xhat);

% known xhat
% xhat = X ;

%% s=[0.5,0.25]
% s1 = 0.5 * ones(N, 1); 
% zhat = A * s1 + Hd ;
% xhat = (zhat' * Y)/(zhat' * zhat); 
% X_sam = zeros(length(Sam), L) ;
% for i = 1 : length(Sam)
%        X_sam(i,:) = abs( Sam(i) - xhat );
% end
% [~,I] = min (X_sam);
% xhat = Sam(I);
%[~ , Xerr] = symerr(X , xhat)


%% BiGAMP estimate xhat

    opt = BiGAMPOpt;
    error_function = @(qval) 20*log10(norm(qval - Z,'fro') / norm(Z,'fro'));
    opt.error_function = error_function;

    [xhat] = BiGAMP(opt, PrioriIn, optIn );

    phase = Sam(1)/xhat(1) ;
    xhat = phase * xhat ;
    X_sam = zeros(length(Sam), L) ;
    for i = 1 : length(Sam)
           X_sam(i,:) = abs( Sam(i) - xhat );
    end
    [~,I] = min (X_sam);
    xhat = Sam(I);
   %  [~ , Xerr] = symerr(X , xhat);
     
%% OMP
%Zhat = (Y * xhat')/L - Hd ;
%[~,shat] = SOMP2(Zhat, A, floor(rho*N) );
%[model,~] = rvmRegVb(A.',Zhat.');
%[shat, t] = cosamp(A,Zhat, ceil(rho*N), 1e-3,50);
%% GAMP

    Y_s = (Y * xhat')/L - Hd ;
    nuw_s = nuw/L ;
    [shat ] = GAMP_Bernoulli( Y_s, A, nuw_s, rho, s ); 
  %  [shat ] = GAMP_Bernoulli( Y_s, A, nuw_s, rho, s ); 
 
%%
% known shat
% shat = s ;
% zhat = A * shat + Hd ;
% xhat = (zhat' * Y)/(zhat' * zhat); 

results.xhat = xhat;
results.shat = shat;

end
