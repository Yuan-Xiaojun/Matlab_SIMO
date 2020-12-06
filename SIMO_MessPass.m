function results = SIMO_MessPass( optIn, PrioriIn )

%%
M = optIn.M;
N = optIn.N;
L = optIn.L;
nit = optIn.nit;
QAM = optIn.QAM;


nuw = PrioriIn.noiseVar ;
nuX = PrioriIn.nuX ;
A = PrioriIn.A;
Hd = PrioriIn.Hd;
Y = PrioriIn.Y;
Sam = PrioriIn.Sam;
State = PrioriIn.State;

% Y = pinv(A)*Y;
% A = eye(M);

%%
A2 = abs(A).^2;
Hd2 = abs(Hd).^2;

stop = false;

%shat = 0.5*ones(M,L,N);
shat = randi([0,1],M, L, N);
% for  n = 1 : N
%     shat(:,:,n) = PrioriIn.s(n);
% end
% svar = zeros(M,L,N);
shat2 = abs(shat).^2;
svar = 0.25*ones(M,L,N);

zhat = zeros(M,L);
zvar = ones(M,L);
pv1 = zeros(M,L);

qvar = ones(M,L,length(Sam));
qhat = zeros(M,L,length(Sam));
Px2y = zeros(M,L,length(Sam));

rvar = ones(M,L,N,length(State));
rhat = zeros(M,L,N,length(State));
Ps2y = zeros(M,L,N,length(State));

% xhat = repmat(PrioriIn.X, M,1);
% xhat2 = abs(xhat).^2;
% xvar = zeros(M,L);

it = 0;    
while ~stop

    it = it + 1;

    % Check for final iteration
    if it >= nit
        stop = true;
    end

    %%  message passing between check node y and variable node x

    %     % computate zhat zvar
    for m = 1 : M
        zhat(m, :) = A(m, :) * reshape(shat(m,:,:), L, N).' + Hd(m);
        zvar(m, :) = A2(m, :) * reshape(svar(m,:,:), L, N)';
    end
    % computate Py2x
    for i = 1 : length(Sam)
        qvar(:,:,i) = zvar.*abs(Sam(i))^2;
        qhat(:,:,i) = Y - zhat.*Sam(i);
    end
    qvar = qvar + nuw;
    Py2x = - (log(qvar) + abs(qhat).^2 ./qvar);      % log  version  log(Py2x)
    Px = sum(Py2x);                                    % posterior of x
    for m = 1 : M
        Px2y(m,:,:) = Px(1,:,:) - Py2x(m,:,:);            % log  version  log(Px2y)
    end
    %because the value of Px2y is too small, 
    %so we multiply a mutual factor to provent the case of "NAN" arise.
    Px2y_min =  max(Px2y,[], 3);
    for i = 1 : length(Sam)
            Px2y(:, :, i) = Px2y(:, :, i) - Px2y_min;
    end
    Px2y = exp(Px2y);     
    % normailzed Px2y
    P_x = zeros(M,L);
    for i = 1 : length(Sam)
        P_x = P_x + Px2y(:,:,i);
    end
    % P_x = max(P_x, eps);
    for i = 1 : length(Sam)
        Px2y(:,:,i) = Px2y(:,:,i)./P_x;
    end
    % Compute the mean and variance  of Px2y  
    xhat = zeros(M,L);
    for i = 1 : length(Sam)
       xhati = Px2y(:,:,i).*(Sam(i).*ones(M,L));
       xhat = xhat + xhati;              
    end
    xvar = zeros(M,L);
    for i = 1 : length(Sam)
       xvari = Px2y(:,:,i).*(abs(Sam(i).*ones(M, L) - xhat).^2);
       xvar = xvar + xvari;              
    end
    xhat2 = abs(xhat).^2;

    %%  message passing between check node y and variable node s

    % computate phat pvar
    for m = 1 : M
        pv1(m, :) = A2(m, :) * reshape(shat2(m,:,:) + svar(m,:,:) , L, N)' + Hd2(m);
    end
    phat = zhat .* xhat;
    pvar = pv1.* xvar + zvar.*xhat2;
    % expand demension to dot product
    YN = repmat(Y, 1,1,N);
    phatN = repmat(phat, 1,1,N);
    pvarN = repmat(pvar, 1,1,N);
    xhatN = repmat(xhat, 1,1,N);
    xvarN = repmat(xvar, 1,1,N);
    xhat2N = repmat(xhat2, 1,1,N);
    AL = reshape(repmat(A,L,1), M,L,N);
    A2L = reshape(repmat(A2,L,1), M,L,N);

    % computate Py2s
    for i = 1 : length(State)
        rvar(:,:,:,i) = nuw + pvarN + A2L.*(abs(State(i))^2 ...
            - shat2 - svar).*xvarN - A2L.*svar.*xhat2N;

        rhat(:,:,:,i) = YN - phatN - AL.*(State(i) - shat).*xhatN;
    end
    Py2s = - (log(rvar) + abs(rhat).^2 ./rvar);      % log  version  log(Py2s)
    Ps = sum(sum(Py2s, 1), 2);                                    % posterior of s
    for m = 1 : M
        for l = 1 : L
            Ps2y(m,l,:,:) = Ps(1,1,:,:) - Py2s(m,l,:,:);             % log  version  log(Ps2y)
        end
    end  

    %Multiply a mutual factor to provent the case of "NAN" arise.
    %Ps2y = Ps2y/(M*L) ;
    Ps2y_min =  max(Ps2y,[], 4);
    for i = 1 : length(State)
            Ps2y(:, :, :, i) = Ps2y(:, :, :, i) - Ps2y_min;
    end
     Ps2y = exp(Ps2y) ; 

    % normailzed Ps2y
    P_s = zeros(M,L,N);
    for i = 1 : length(State)
        P_s = P_s + Ps2y(:,:,:,i);
    end
    %P_s = max(P_s, eps);
    for i = 1 : length(State)
        Ps2y(:,:,:,i) = Ps2y(:,:,:,i)./P_s;
    end
    % Compute the mean and variance  of Ps2y  
    shat = zeros(M,L,N);   
    for i = 1 : length(State)
       shati = Ps2y(:,:,:,i).*(State(i).*ones(M,L,N));
       shat = shat + shati;              
    end
    svar = zeros(M,L,N);
    for i = 1 : length(State)
       svari = Ps2y(:,:,:,i).*(abs(State(i).*ones(M,L,N) - shat).^2);
       svar = svar + svari;              
    end 
    shat2 = abs(shat).^2;

end
results.Px = reshape(Px(1,:,:), L, length(Sam))';
results.Ps = reshape(Ps(1,1,:,:), N, length(State))';  

end
