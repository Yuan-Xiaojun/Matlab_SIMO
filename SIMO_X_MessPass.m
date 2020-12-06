function [P_x2y1, P_x] = SIMO_X_MessPass( Zhat,Zvar, optIn,  PrioriIn )

%%
M = optIn.M;
N = optIn.N;
L = optIn.L;
nit = optIn.nit;

nuw = PrioriIn.noiseVar ;
Hd = PrioriIn.Hd;
Y = PrioriIn.Y;
Sam = PrioriIn.Sam;
C = length(Sam);
X = PrioriIn.X;

%%  Initailization
Sam2 = abs(Sam).^2 ;

stop = false;

P_x2y = ones(M, L, C)./C;     
xvar_ml = ones(M, L, C);
xhat_ml = ones(M, L, C);

%%
it = 0;    
while ~stop

    it = it + 1;

    % Check for final iteration
    if it >= nit.innerx
        stop = true;
    end
%%

     % message from y_ml to x_l
     for i = 1 : C
         xvar_ml(:,:,i) = nuw + Zvar.ml * Sam2(i) ;
         xhat_ml(:,:,i) = Y - (Zhat.ml + repmat(Hd, 1, L))* Sam(i) ;
     end
     P_y2x = - (log(xvar_ml) + abs(xhat_ml).^2 ./xvar_ml);      % log  version  log(P_y2x)
     % message from x_l to y_ml
     P_x = sum(P_y2x, 1);                                    % posterior of x
     for m = 1 : M
         P_x2y(m,:,:) = P_x(1,:,:) - P_y2x(m,:,:);            % log  version  log(Px2y)
     end
     % Multiply a mutual factor to provent the case of "NAN" arise.
     Px2y_max =  max(P_x2y,[], 3);   
     for i = 1 : length(Sam)
         P_x2y(:, :, i) = P_x2y(:,  :, i) - Px2y_max;
     end
     P_x2y = exp(P_x2y) ; 
     % Normolization
     P_sum = zeros(M,L);
     for i = 1 : length(Sam)
         P_sum = P_sum + P_x2y(:,:,i);
     end
     % P_x = max(P_x, eps);
     for i = 1 : length(Sam)
         P_x2y(:,:,i) = P_x2y(:,:,i)./P_sum;
     end 
    
end 
    [~,b] = max(P_x2y(1,1,:));   %  eliminate  ambiguity
    P_x2y1 = zeros( size(P_x2y));
    for i = 1 : C
        col = mod((i + b - 1),  4);
        if col == 0
            col = 4;
        end
        P_x2y1(:,:,i) = P_x2y(:,:, col );
    end
    
%  xhat = zeros(M,N);   
%  for i = 1 : length(Sam)
%      xhati = P_x2y1(:,:,i).*(Sam(i).*ones(M,N));
%      xhat = xhat + xhati;              
%  end
%  xvar = zeros(M,N);
%  for i = 1 : length(Sam)
%     xvari = P_x2y1(:,:,i).*(abs(Sam(i).*ones(M,N) - xhat).^2);
%     xvar = xvar + xvari;              
%  end 
end
