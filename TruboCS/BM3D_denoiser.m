%denoiser using BM3D
%you should install bm3d first
function [F,F_div] = BM3D_denoiser(y,v)
[m,n]=size(y);
% if m~=n
%     if floor(sqrt(m*n))-sqrt(m*n)~=0
%     else
%         n=sqrt(m*n);
%         k=reshape(y,n,n);
%     end
% else
%     k=y;
% end

n=sqrt(m*n);
% y=reshape(y,n,n);
denoiser='fast-BM3D';
%denoiser='BM3D';
denoi=@(noisy,sigma_hat) denoise(noisy,sigma_hat,n,n,denoiser);

% denoised=denoi(y,v);
% eta=randn(1,n^2);
% epsilon=max(y)/1000+eps;
% div=eta*((denoi(y+epsilon*eta',v)-denoised)/epsilon)/n^2;

% [~,y_hat]=BM3D(1,k,v,'np',0);
% xhat = obj.fxnDenoise(rhat,rvar);
denoised=denoi(y,v);
div=0;
K=1;
for ii=1:K
    %epsilon = 0.05*min(v, mean(abs(y),1)) + eps;
    epsilon=max(y)/100+eps;
    %eta = sign(randn(1,n^2)); % random direction
    eta=randn(1,n^2);
    rhat_perturbed = y + epsilon*eta';
    xhat_perturbed=denoi(rhat_perturbed,v);
    divp=eta*(xhat_perturbed-denoised)/epsilon;
    div=div+divp;
end
div=div/K/n^2;


% de_type='lc';%lc,np
% [~,y_hat]=BM3D(1,k,v,de_type,0); %fast BM3D
% y_hat=y_hat*255;
% eta=randn(n,n);
% epsilon=max(y)/1000+eps;
% [~,y_hat1]=BM3D(1,k+epsilon*eta,v,de_type,0);
% y_hat1=y_hat1*255;
% F_div=sum(sum(eta.*((y_hat1-y_hat)/epsilon)))/n^2;
% F=reshape(y_hat,1,n^2);

F=reshape(denoised,1,n^2);
F_div=div;
end