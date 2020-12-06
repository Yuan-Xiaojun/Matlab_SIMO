function [xhat] = GAMP_S(Y, A, nuw, rho, s  )

%%
[M, T] = size(Y);
[~, K] = size(A);
State = [0, 1] ;
bata_shat = 0.1;
bata_x = 0.2;

%mean(A(:))
%%
%gX = CAXwgnEstimIn(0,1);
gOut = CAwgnEstimOut(Y, nuw);
%gX =  AwgnEstimOut(0.5, 0.25);

xvar = rho*(1-rho)*ones(K, T);
xhat= rho*ones(K, T);
% xhat = ones(K,1);
% p = randperm(K) ;
% xhat(p(1: floor(K * (1 - rho)))) = 0;

shat=zeros(size(Y)); % used for the tailor series  expansion
Ahat2 = abs(A).^2;
xvarMin=0;
zvarToPvarMax=0.99;
N_amp=100;
amp_diff_stop = 1e-10 ;

svar = 0 ;
rvar = 0 ;
pvar = 0 ;
rhat = 0 ;


for it_gamp=1:N_amp
     zvar = Ahat2*xvar;
     
  %   pvar_last = pvar ;
     pvar = zvar;
  %   pvar = (1 - bata_shat) * pvar_last + bata_shat * pvar ;  % damping
     
     zhat = A * xhat;
     phat = zhat- shat.*zvar;
     [zhat0,zvar0] = gOut.estim(phat,pvar);
     pvarInv = 1 ./ pvar;
     shat_last = shat ;
     shat = pvarInv.*(zhat0-phat);
     shat = (1 - bata_shat) * shat_last + bata_shat * shat ;  % damping
     
     svar_last = svar ;
     svar = pvarInv.*(1-min(zvar0./pvar,zvarToPvarMax));
     svar = (1 - bata_shat) * svar_last + bata_shat * svar ;  % damping
     
     rvar_last = rvar ;
     rvar = 1./(Ahat2.'*svar);
     rvar = (1 - bata_shat) * rvar_last + bata_shat * rvar ;  % damping
     
     rhat_last = rhat ;
     rhat = xhat + rvar.*(A'*shat);
     rhat = (1 - bata_shat) * rhat_last + bata_shat * rhat ;  % damping
     
     rvar = max(rvar, xvarMin);
     xhat_old = xhat ;
     [xhat,xvar ] = Xestim(rhat, rvar );
   %  [xhat,xvar, ~ ] = gX.estim( rhat, rvar ) ;
     xhat = (1 - bata_x) * xhat_old + bata_x * xhat ;  % damping 
     if norm(xhat - xhat_old)^2/numel(xhat) < amp_diff_stop
         break; 
     end
end

    X_sam = zeros(K, T, length(State)) ;
    for i = 1 : length(State)
        X_sam(:,:, i) = abs( State(i) - xhat );
    end
    [~, I] = min (X_sam, [], 3);
    xhat = State(I);
    
%     xhat = xhat(1:end-2); 
%     xvar = xvar(1:end-2); 

end

