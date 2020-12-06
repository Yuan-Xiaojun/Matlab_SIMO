% --------------------------
clear;
clc;
%--------------------------
rou=0.1;
N=3000;
delta=0.5;
M=delta*N;
M=fix(M);
sigma=sqrt(10^(-50/10));
% sigma=0;
%=========distribution of x=========
pattern=random('Binomial',1,rou,N,1);
x = 1/sqrt(rou)*randn(N,1);
x = x.* pattern;
x = x/norm(x)*sqrt(N);
permuation=randperm(N);
index_selection=permuation(1:M);
p=2*(rand(N,1)>0.5)-1;
% p2=2*(rand(N,1)>0.5)-1;
% xtmp=1/sqrt(delta)*dct((p.*idct(x))); % xtmp= F*x;
xtmp = 1/sqrt(delta)*dct((p.*idct(x)));
y=xtmp(index_selection)+sigma*randn(M,1); %y = S*xtmp + w (S*F = AS)

iter=30;
K=3;
epsilon=1e-6;
errx = @(z) norm(z-x)^2/norm(x);
[errx_out,x_hat] = dturbocs(y,p,M,N,K,index_selection,sigma,errx);
semilogy(errx_out)





