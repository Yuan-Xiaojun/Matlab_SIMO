function [ s_est ] = GAMP_Bernoulli( Y, A, nuw, rho, s )

% This is a message passing algorithm for compressed sensing problems with discrete variables (0,1). 
% The algorithm doesn't use these approximations in GAMP when combining the messages from check nodes to variable nodes,
% but calculate these messages directly. So it slihtly outperforms the traditional GAMP in most cases. 
% Paramater setting: Y = Z + W = As + W ; 
% Y: observed signal;
% A: sensing matrix;
% s: true value;
% W: noise~CN(0, nuw);
%rho: the probability of s taken 1 ;
% shat.m2n and svar.m2n  are the mean and variance of s passed from check nodes to variable nodes
% shat.m2n and svar.m2n are M*N*length(State) matrixes; the (m,n,i)th entry is the message from y_m to s_n when s_n taken as its i-th value(i.e., 0 or 1).
% shat.n2m and svar.n2m are the mean and variance of s passed from variable nodes to check nodes.
% shat.m2n and svar.m2n are M*N matrixes. the (m,n)th entry is the message from s_n to y_m. 

%%
if nargin == 0
    clc
    %Problem dimensions
    optIn.M = 20; %number of observation signals
    optIn.N = 20; %number of detect signals
    rho = 0.5; 
    
    A = sqrt(0.5)*(randn(optIn.M, optIn.N) + 1i*randn(optIn.M, optIn.N)) ;
    s = randi([0,1], optIn.N, 1);
    nuw = 1e-4 ;
    Y = A * s + nuw ;   
end

[M, N] = size(A) ;
State = [0, 1];
Prior_s = rho/(1-rho); % P_1/P_0
% Initialization
A2 = abs(A).^2 ;
shat_m2n = zeros(M,N,length(State)) ;
P_s2y = zeros( M, N, length(State));
shat_n2m = ones(M,N);
for m = 1 : M
    p = randperm(N) ;
    shat_n2m(m , p(1: floor(N * (1 - rho)))) = 0;
end
svar_n2m = (rho - rho^2)*ones(M, N);
Zvar = sum(A2 .* svar_n2m , 2) ;
Zhat = sum(A .* shat_n2m , 2);

stop = false;
it = 1 ;
nit = 20;
amp_stop = 1e-8 ;
%%
while ~stop
%%
 % message from check node to s_n
 % repmat(Y-Zhat,1,N)+A.*shat.n2m is the message from y_m to s_n after remove the exterior message from s_n to y_m.
 % [(Y-Zhat,1,N)+A.*shat.n2m]_{m,n} = A_{m,i}*s_i + w + As_{tilde};(s_{tilde}~CN(0,svar.n2m))
 % the mean of A_{m,i}*s_i is [(Y-Zhat,1,N)+A.*shat.n2m]_{m,n}. 
     for i = 1 : length(State)
         shat_m2n(:,:,i) = repmat(Y - Zhat , 1, N) - A .*( State(i) - shat_n2m) ;
     end   
     svar_m2n = repmat( (repmat( nuw + Zvar, 1, N) - A2 .* svar_n2m ), 1, 1, length(State));
     svar_m2n = max(svar_m2n, 1e-100); 
     P_y2s = -  abs(shat_m2n).^2 ./svar_m2n ;      % log  version  log(P_dletas) 
     % message from s_n to delta
     P_s = sum(P_y2s ,1);   % posterior of s   
     for m = 1 : M
         P_s2y(m,:,:) = P_s(1,:,:) - P_y2s(m,:,:); % remove exterior message
     end 
     % avoid data overflow
%      Max_Ps2y = max(P_s2y,[], 3) ;
%      for i = 1 : N  
%          P_s2y( Max_Ps2y(:,i) < -50, i, :) = -1 ;
%      end        
     P1_P0 = P_s2y(:,:,2) - P_s2y(:,:,1) ;
     P1_P0 = max(min(P1_P0, 50 ), -50);
     P10_s2y = exp( P1_P0 ) ; % P1_s2y/P0_s2y
     %  multily prior
     P10_s2y = P10_s2y * Prior_s ;
     % Normalization 
     P_s2y (:,:,2) = P10_s2y ./(1 + P10_s2y); 
     P_s2y(isnan(P_s2y)) = 1;  % eps/eps=1; inf/inf =NaN
     P_s2y (:,:,1) = 1 - P_s2y (:,:,2);
     % Compute the mean and variance of message from s_n to y_m  
     shat_old = shat_n2m ;
     shat_n2m = P_s2y (:,:,2) ; %the expectation of s is the probability of taking 1.
     svar_n2m = P_s2y(:,:,1).* P_s2y (:,:,2) ; %the expectation of s is rho*(1-rho).
     % Compute the mean and variance  of z = As ; 
     Zvar = sum(A2 .* svar_n2m , 2) ;
     Zhat = sum(A .* shat_n2m , 2);   
    
     % Check for final iteration
    it = it + 1;
%     if (it >= nit)||(norm(shat_n2m - shat_old, 'fro')^2/M/N <= amp_stop)
%         stop = true;
%     end   
    if (it >= nit)
        stop = true;
    end  


%      Ps = reshape(P_s(1,:,:), N, length(State)); 
%      [~,I] = max( Ps,[],2 ) ;
%      s_est = State( I )' ;
%      [~ , Serr] = biterr(s , s_est);
%      Serr
%      norm(Y-A*s)^2
end
log_Ps = reshape(P_s(1,:,:), N, length(State));  
% Ps10 = max(min( ( log_Ps(:,2) - log_Ps(:,1)) , 50 ), -50);
% Ps10 = exp( Ps10) ; % P_1/P_0
% shat = Ps10./(1 + Ps10); % the expectation of s is the probability of taking 1.
% shat(isnan(shat)) = 1;
% svar = shat.*(1-shat) ;

 [~,I] = max( log_Ps,[],2 ) ;
 s_est = State( I )' ;
%  [~ , Serr] = biterr(s , s_est);
%  Serr
%  norm(Y-A*shat)^2/numel(Y)

end







