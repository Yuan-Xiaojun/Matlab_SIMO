function [shat, svar, P_s] = SIMO_S_MessPass( Zhat, Zvar, shat, svar, optIn,  PrioriIn )

if nargin == 0
    clc

    %Problem dimensions
    optIn.M = 30; %coherence time
    optIn.N = 20; %transmit antenna
    optIn.L = 20; %received antenna why set in this way?
    
    PrioriIn.Sp_s = 1;
    PrioriIn.State = [0, 1];
    
    PrioriIn.A = sqrt(0.5)*(randn(optIn.M, optIn.N) + 1i*randn(optIn.M, optIn.N)) ;
    PrioriIn.s = randi([0,1], optIn.N, 1);
    Zvar.mdelta = 1e-4 ;
    Zhat.mdelta = PrioriIn.A * PrioriIn.s + Zvar.mdelta ;   
    optIn.nit.inner = 10 ;
end
M = optIn.M;
N = optIn.N;
nit = optIn.nit;

A = PrioriIn.A;
State = PrioriIn.State;
Sp_s = PrioriIn.Sp_s;
s = PrioriIn.s;
rho = PrioriIn.rho ;
A2 = abs(A).^2 ;
shat.mn = zeros(M,N,length(State)) ;
P_s2dleta = zeros( M, N, length(State));

shat.nm = ones(M,N);
for m = 1 : M
    p = randperm(N) ;
    shat.nm(m , p(1: floor(N * (1 - rho)))) = 0;
end
svar.nm = (rho - rho^2)*ones(M, N);

Zvar.deltam = sum(A2 .* svar.nm , 2) ;
Zhat.deltam = sum(A .* shat.nm , 2);

stop = false;
it = 0 ;
%%
while ~stop
    
    it = it + 1;

    % Check for final iteration
    if it >= nit.inners
        stop = true;
    end
%%
 % message from delta to s_n
     for i = 1 : length(State)
         shat.mn(:,:,i) = repmat(Zhat.mdelta - Zhat.deltam , 1, N) - A .*( State(i) - shat.nm) ;
     end    
     svar.mn = repmat( (repmat(Zvar.mdelta + Zvar.deltam , 1, N) - A2 .* svar.nm ), 1, 1, length(State));
     svar.mn = max(svar.mn, 1e-200); 
     P_dleta2s = - (log(svar.mn) + abs(shat.mn).^2 ./svar.mn);      % log  version  log(P_dletas) 
     % message from s_n to delta
     P_s = sum(P_dleta2s ,1);                                % posterior of s   
     for m = 1 : M
         P_s2dleta(m,:,:) = P_s(1,:,:) - P_dleta2s(m,:,:);             % log  version  log(P_sdleta)
     end 
     P10_s2dleta = min(max(exp( P_s2dleta(:,:,2) - P_s2dleta(:,:,1) ), eps), inf);
     %  multily prior
     P10_s2dleta = P10_s2dleta * Sp_s ;
     % Normolization 
     P_s2dleta (:,:,2) = P10_s2dleta ./(1 + P10_s2dleta);
     P_s2dleta(isnan(P_s2dleta)) = 1;
     P_s2dleta (:,:,1) = 1 - P_s2dleta (:,:,2);
     % Compute the mean and variance  of P_sdleta  
     shat.nm = zeros(M,N);   
     for i = 1 : length(State)
         shati = P_s2dleta(:,:,i).*(State(i).*ones(M,N));
         shat.nm = shat.nm + shati;              
     end
     svar.nm = zeros(M,N);
     for i = 1 : length(State)
        svari = P_s2dleta(:,:,i).*(abs(State(i).*ones(M,N) - shat.nm).^2);
        svar.nm = svar.nm + svari;              
     end 
     % Compute the mean and variance  of delta to m  
     Zvar.deltam = sum(A2 .* svar.nm , 2) ;
     Zhat.deltam = sum(A .* shat.nm , 2);
     
%      Ps = reshape(P_s(1,:,:), N, length(State))'; 
%      [~,I] = max( Ps ) ;
%     sh = State( I ) ;
%     cor_s = length(find((s - sh' ) ==0));
%     Serr = (N - cor_s)/N
%      
end
%      Ps = reshape(P_s(1,:,:), N, length(State))'; 
%      [~,I] = max( Ps ) ;
%     sh = State( I ) ;
%     cor_s = length(find((s - sh' ) ==0));
%     Serr = (N - cor_s)/N 
end