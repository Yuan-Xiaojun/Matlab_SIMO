function results = SIMO_algorithm(optIn)


%%
M = optIn.M;
N = optIn.N;
L = optIn.L;
SNR = optIn.SNR;
QAM = optIn.QAM;
State = [0,1];
rho = optIn.rho;
Sp_s = rho/(1-rho) ;    % P1/P0
beta = 0.5;   % fadding factor
%% produce System model 

if QAM == 4
    Sam = [1+1i -1+1i -1-1i 1-1i]./sqrt(2);    
elseif QAM == 16
    Sam = [1+1i -1+1i -1-1i 1-1i 3+1i -1+3i -3-1i 1-3i 3+3i -3+3i -3-3i 3-3i 1+3i -3+1i -1-3i 3-1i]./sqrt(10);
end

G = sqrt(0.5)*(randn(M,N) + 1i*randn(M,N));
Hd = sqrt(0.5)*(randn(M, 1) + 1i*randn(M, 1));
Hr = sqrt(0.5)*(randn(N, 1) + 1i*randn(N, 1));
Dh = diag(Hr);
theta = CVX_Optimal( G, Hd, Dh, rho, N);
%theta = randn(N,1) + 1i*randn(N,1);
%theta = theta ./(abs(theta));
A = beta * G * diag(theta) * Dh ; 
% Produce RIS data
s = ones(N,1);
p = randperm(N) ;
s(p(1: floor(N * (1 - rho)))) = 0;
% Produce user signal
B = randi([0,1], sqrt(QAM), L);
B(:,1) = 1;
% X = bi2de(B','left-msb');
% X = qammod(X, QAM).';
% X = X./abs(X(1)) ;
X = Constell_Modulate(B, QAM);
nuz = 1 + N * beta* rho * 0.5 ;
nuw = 10^(-SNR/10) ;
Z = (A*s + Hd) * X ;
Y = Z + sqrt(nuw/2)*(randn(size(Z))+1i*randn(size(Z))) ;

%% Initialize results as empty
results = [];
PrioriIn.noiseVar = nuw;
PrioriIn.ZVar = nuz;
PrioriIn.nuX = 1;
PrioriIn.A = A;
PrioriIn.X = X;
PrioriIn.Hd = Hd;
PrioriIn.Y = Y;
PrioriIn.Sam = Sam;
PrioriIn.QAM = QAM;
PrioriIn.State = State ;
PrioriIn.Sp_s = Sp_s ;
PrioriIn.rho = rho ;

PrioriIn.s = s;
PrioriIn.X = X;
%% 

   % [results1] = SIMO_MessPass( optIn, PrioriIn ); 
    %[results1] = SIMO_test( optIn, PrioriIn );
    [results1] = SIMO_SVD( optIn, PrioriIn );

    xhat = results1.xhat ;
    shat = results1.shat ;
  

%%
     Bhat = Constell_Mapping(xhat, QAM, Sam) ;    
%     Bhat = qamdemod(xhat, QAM);
%     Bhat = de2bi(Bhat,'left-msb');   
     [~ , Xerr] = biterr(B , Bhat);
     [~ , Serr] = biterr(s , shat );
     
     results.Xerr = Xerr ;
     results.Serr = Serr ;


    
end
%%







