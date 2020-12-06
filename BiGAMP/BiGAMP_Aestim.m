function [Ahat,Avar,val] = BiGAMP_Aestim(rhat, rvar, nuz)


uhat0 = 0;
uvar0 = nuz; 

% Compute posterior mean and variance
gain = uvar0./(uvar0+rvar);
Ahat = gain.*(rhat-uhat0)+uhat0;
Avar = gain.*rvar;

u1 = sum(Ahat(:))/numel(Ahat);
v1 = sum(sum( abs(Ahat - u1*ones(size(Ahat))).^2 ))/numel(Ahat)  +sum(Avar(:))/numel(Ahat);
val =  (log(v1./uvar0(1,1)) + (1-v1./uvar0(1,1)) ...
    - abs(u1-uhat0(1,1)).^2./uvar0(1,1) );