% W.J.Yan, 3.2019

clc;
clear;
tic
% randn('state',1);
% rand('state',1);
% randi('state',1);
setup_DL 

N_Rec = [32];   % receiver antenna
N_IRS = [16];   % IRS antenna
N_L = 100;       % coherencr time
SNR = [-20 : 1 : -10]; 
SimLen = 5000;   % average simulation time

rho = 0.9 ;  % mean_s
nit.outer = 200 ;   % limit iterations for SIMO_Z_MessPass 
nit.inners = 10;   % limit iterations for SIMO_S_MessPass  
nit.innerx = 3 ;
QAM = 4;   % modulate style,  only for: 4 : QPSK;  16 : 16QAM.


for iR = 1 : length(N_Rec)
    for iS = 1 : length(N_IRS)
        
        optIn.M = N_Rec(iR);    
        optIn.N = N_IRS(iS);
        optIn.L = N_L;      
        optIn.QAM = QAM;
        optIn.rho = rho;
        optIn.nit = nit; %limit iterations
        
        X_BER = zeros(1, length(SNR));    % average of bit error rate of X 
        S_BER = zeros(1, length(SNR));  % average of bit error rate of S
        for sim = 1 : SimLen
            
            tX_BER = zeros(1, length(SNR));    % bit error rate of X in each iteration 
            tS_BER = zeros(1, length(SNR));  % bit error rate of S in each iteration   
%%            
            for iN = 1 : length(SNR)  
                
                optIn.SNR = SNR(iN);
                
                results = SIMO_algorithm(optIn);   %  subfunction
                
                tX_BER(iN) = results.Xerr;
                tS_BER(iN) = results.Serr;
  
            end
%%            
            X_BER = X_BER*(sim-1)/sim + tX_BER/sim;
            S_BER = S_BER*(sim-1)/sim + tS_BER/sim;
             sim
%             fprintf('sim = %3d, M = %3d, N = %3d, L = %3d,  QAM = %2d, rho = %2d\n', sim, optIn.M, optIn.N, optIn.L, optIn.QAM, rho );
%                     
%             fprintf('SNR =');
%             for iN = 1:length(SNR)
%                 fprintf('%2d  ', SNR(iN));
%             end
%             fprintf('\n');
%             
%             fprintf('tX_BER =');
%             for iN = 1:length(SNR)
%                 fprintf('%3.6e  ', tX_BER(iN));
%             end
%             fprintf('\n');
%             
%             fprintf('X_BER =');
%             for iN = 1:length(SNR)
%                 fprintf('%3.6e  ', X_BER(iN));
%             end
%             fprintf('\n');
%             
%             fprintf('tS_BER =');
%             for iN = 1:length(SNR)
%                 fprintf('%3.6e  ', tS_BER(iN));
%             end
%             fprintf('\n');
%             
%             fprintf('S_BER =');
%             for iN = 1:length(SNR)
%                 fprintf('%3.6e  ', S_BER(iN));
%             end
%             fprintf('\n');
            
        end
        
         fid = fopen('SIMO_result.txt','a+');
         fprintf(fid,'\n');
         fprintf(fid,'SimLen = %3d, M = %3d, N = %3d, L = %3d,  QAM = %2d\n', sim, optIn.M, optIn.N, optIn.L, optIn.QAM );
            fprintf(fid, 'SNR =');
            for iN = 1:length(SNR)
                fprintf(fid, '%2d  ', SNR(iN));
            end
            fprintf(fid, '\n');

            fprintf(fid, 'X_BER =');
            for iN = 1:length(SNR)
                fprintf(fid, '%3.6e  ', X_BER(iN));
            end
            fprintf(fid, '\n');

            fprintf(fid, 'S_BER =');
            for iN = 1:length(SNR)
                fprintf(fid, '%3.6e  ', S_BER(iN));
            end
            fprintf(fid, '\n');                 
    end
end
save 
toc
%%






