function results = SIMO_test( optIn, PrioriIn )

%%
M = optIn.M;
N = optIn.N;
L = optIn.L;
nit = optIn.nit;

nuw = PrioriIn.noiseVar ;
A = PrioriIn.A;
Hd = PrioriIn.Hd;
Y = PrioriIn.Y;
Sam = PrioriIn.Sam;
C = length(Sam);
State = PrioriIn.State;
s = PrioriIn.s;
X = PrioriIn.X;
Z = A*s ;
% Step = 0.1;
% StepInc = 1.1;
% StepDec = 0.5;
% StepMin = 0.01;
% StepMax = 0.99;
% DampOpt = 0;
%%  Initailization
A2 = abs(A).^2;
Sam2 = abs(Sam).^2 ;

stop = false;

P_x2y = ones(M, L, C)./C;     % Initailize p_{l,m}^i
% P_x2y = zeros(M, L, C);
% for l = 1 : L
%     [~,c] = min(X(l) - Sam);
%     P_x2y(:,l,c) = 1 ;
% end

% shat.nm = 0.5*ones(M, N);
% svar.nm = ones(M, N);
shat.nm = randi([0,1], M, N);
svar.nm = 0.25*ones(M, N);
% shat.nm = repmat(s', M, 1);
% svar.nm = zeros(M, N);

% Initailize Z_mean, Z_var  keep static during the algorithm.
Z_mean = zeros(M,L,C);
Z_var = nuw*ones(M,L,C);
for i = 1 : C
    Z_mean(:,:,i) = Y/Sam(i) ;
    Z_var(:,:,i) = Z_var(:,:,i)/Sam2(i) ;
end
for m = 1 : M
    Z_mean(m,:,:) = Z_mean(m,:,:) - Hd(m) ;
end
%  Initailize Zvar.deltam, Zhat.deltam.
Zvar.deltam = sum(A2 .* svar.nm , 2);
Zhat.deltam = sum(A .* shat.nm , 2);


%%
it = 0;    
while ~stop

    it = it + 1;

    % Check for final iteration
    if it >= nit.outer
        stop = true;
    end
%%
     % message from y_ml to Z_m
     for i = 1 : C
         P_x2y(:,:,i) = P_x2y(:,:,i)/Sam2(i) ;
     end
     Pxs = sum(P_x2y, 3);    
     for i = 1 : C
         P_x2y(:,:,i) = P_x2y(:,:,i)./Pxs ;
     end
     
     Zhat.lm = reshape(sum(P_x2y.*Z_mean, 3), M, L);   
     Zvar.lm = abs(Z_mean - repmat(Zhat.lm, 1, 1, C)).^2  + Z_var;
     Zvar.lm = reshape(sum(P_x2y.*Zvar.lm, 3), M, L);
     Zvar.lm1 = 1./Zvar.lm ;
     Zmse.lm = abs(norm(Z - Zhat.lm,'fro'))^2 /abs(norm(Z,'fro'))^2 ;
     % message from z_m to delta
     Zvar.mdelta1 = reshape(sum(Zvar.lm1, 2), M, 1) ;
     Zhat.mdelta = reshape(sum(Zhat.lm .* Zvar.lm1, 2), M, 1)./Zvar.mdelta1 ;
     Zvar.mdelta = 1./Zvar.mdelta1 ;
     Zmse.mdelta = abs(norm(Z - Zhat.mdelta,'fro'))^2 /abs(norm(Z,'fro'))^2 ;
%%   
     
%      if (it > 1)
%          Damp = abs(norm(Y - (Zhat.deltam + Hd )* (sum(xhat, 1)/M) ,'fro'))^2 ; 
%          + sum(Zvar.deltam(:)) + sum(xvar(:))/M
%          if ( Damp > DampOpt ) 
%              Step = Step * StepDec ;
%          else
%              Step = Step * StepInc ;
%          end 
%          Step = min( max(Step, StepMin ), StepMax);
%          DampOpt =  Damp ;
%          if it == 2
%              Step = 0;
%          end
%          if it == 3
%              Step = 1;
%          end
%          if Damp < nuw*M*L 
%              stop = true;
%          end       
%          Zvar.deltam = Zvar.deltam * Step + Zvar.deltamOpt * (1 - Step) ;
%          Zhat.deltam = Zhat.deltam * Step + Zhat.deltamOpt * (1 - Step) ; 
%      end
 %%  GAMP module to estimate s   
      
      [shat, svar, P_s] = SIMO_S_MessPass( Zhat, Zvar, shat, svar, optIn,  PrioriIn );  
                    
%%     
    % message from delta to z_m
     Zvar.deltamOpt = Zvar.deltam ;
     Zhat.deltamOpt = Zhat.deltam ;
     
     Zvar.deltam = sum(A2 .* svar.nm , 2);     
     Zhat.deltam = sum(A .* shat.nm , 2);
     Zvar.deltam1 = min(1./Zvar.deltam, inf) ;
     Zmse.deltam = abs(norm(Z - Zhat.deltam,'fro'))^2 /abs(norm(Z,'fro'))^2 ;     
     % message from z_m to y_ml  
     Zhat.lm1 = Zhat.lm .* Zvar.lm1 ;
     for l = 1 : L
         Zvar.ml1(:,l) = Zvar.mdelta1 - Zvar.lm1(:,l) + Zvar.deltam1 ;
         Zhat.ml(:,l) = ( Zhat.mdelta .* Zvar.mdelta1 - Zhat.lm1(:,l) + Zhat.deltam.*Zvar.deltam1 )./Zvar.ml1(:,l) ;
     end
     Zvar.ml = max(1./Zvar.ml1, 1e-200) ;
     Zmse.ml = abs(norm(Z - Zhat.ml,'fro'))^2 /abs(norm(Z,'fro'))^2 ;
%%  GAMP module to estimate X  

     [P_x2y,  P_x] = SIMO_X_MessPass( Zhat,Zvar, optIn,  PrioriIn ) ;
     P_x = reshape(P_x(1,:,:), L, length(Sam))'; 
     [~,I] = max( P_x ) ;
     xhat = Sam( I ) ;
     [~ , Xerr] = symerr(X , xhat)
end

% P_x = reshape(P_x(1,:,:), L, length(Sam))'; 
% [~,I] = max( P_x ) ;
% xhat = Sam( I ) ;

results.xhat = xhat;
results.Ps = reshape(P_s(1,:,:), N, length(State))'; 

%     [~,I] = max( results.Px ) ;
%     xhat1 = Sam( I ) ;
%     cor_x = length(find((X-xhat1)==0));
%     Xerr = ( L - cor_x)/( L) ;
%     
%      [~,I] = max( results.Ps ) ;
%     sh = State( I ) ;
%     cor_s = length(find((s - sh' ) ==0));
%     Serr = (N - cor_s)/N    

end
