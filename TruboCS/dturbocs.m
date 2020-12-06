% OAMP_KERNEL_DCT FUNCTION (!!!fast computation version!!!)
% According to algorithm OAMP
% this function should be applied to Orthogonal matrix A for we can use
% fast algorithm. We assume the measurement matrix have been row normalized
% to one.
% we should know the row choosen by measurement matrix.
function [errx_out,x_hat]=dturbocs(y,p,M,N,K,index,sigma,errx)
% input signal: y
% measurement matrix: A
% w=At*(A*At)^(-1);
% a=(n^2)/trace(w*A);
% AG=a*w;
% inverse matrix: AG;
% trace of At*A: TRA
%-------------parameters----------------
% dimension of measurement matrix, measurement rate delta
delta=M/N;
% iterative times
iter=50;
% initial x_hat and z_hat
x_hat=zeros(N,1);
z_hat=y;
z_hatpri=zeros(M,1);
x_p=zeros(N,1);
% parameter to ..
epsilon=1e-6;
% noise level, in the case of noiseless
%sigma=0;
%-----------iterative algorithm------------------
for ii=1:iter
    % noise variance of z_hat
    v=(norm(z_hat)^2-M*sigma^2)/N;
    % LE step of OAMP
    ztmp=zeros(N,1);
    ztmp(index)=z_hat;
    r_hat=x_hat+1/sqrt(delta)*dct((p.*idct(ztmp)));
    %r_hat=x_hat+A'*z_hat
    % estimated noise of LE
    c=(1/delta)*v+N/M*sigma^2;
    %c=v;
    % Kernel function to promote sparsity
    if K==3
        [F,F_div]=Kernel_lin_1(r_hat,c); 
    elseif K==2
        [F,F_div]=Kernel_lin_3(r_hat,c);
    elseif K==1
        %--BM3D denoiser used in original image denoise--
        [F,F_div]=BM3D_denoiser(r_hat,sqrt(c));
    end
    F_kernel=zeros(K,N);
    for jj=1:K
        F_kernel(jj,:)=F(jj,:)-F_div(jj)*reshape(r_hat,1,N);
    end
    %---estimate of x by sure minimization----
    C_ff=F_kernel*F_kernel'/N+epsilon*eye(K);
    C_fy=conj(F_kernel)*r_hat/N;
    C_cal=mldivide(C_ff,C_fy);
    % estimated signal
    x_hat=F_kernel.'*C_cal;
    %--------residual of estimated signal------
    xtmp=1/sqrt(delta)*dct((p.*idct(x_hat)));
    z_hat=y-xtmp(index);
    z_hatpri=z_hat;
    
    x_p=x_hat;
    errx_out(ii) = errx(x_hat);
end
% it
%-----------sure let output step----------
P_ff=F*F'/N;
P_ff=P_ff+epsilon*eye(K);
P_fx=1/N*F*r_hat-c*F_div';
a=mldivide(P_ff,P_fx);
x_hat=F'*a;
%--------------------
end