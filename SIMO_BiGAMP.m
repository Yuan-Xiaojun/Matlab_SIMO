function results = SIMO_BiGAMP( optIn, PrioriIn )

%%
M = optIn.M;
N = optIn.N;
L = optIn.L;
nit = optIn.nit;

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
     [~ , Xerr] = symerr(X , xhat);


%% GAMP

    Zhat.mdelta = (Y * xhat')/L - Hd ;
    Zvar.mdelta = nuw/L ;
    [~, ~, P_s ] = SIMO_S_MessPass( Zhat, Zvar, shat, svar, optIn,  PrioriIn ); 
    
%end
%%
% known shat
% shat = s ;
% zhat = A * shat + Hd ;
% xhat = (zhat' * Y)/(zhat' * zhat); 

results.xhat = xhat;
%results.shat = shat;
results.Ps = reshape(P_s(1,:,:), N, length(State))'; 

end
