function [xhat,xvar,vall] = BiGAMP_Xestim(qhat, qvar, Sam)

[U,L]=size(qhat);
xhat=zeros(U,L);
xvar=zeros(U,L);
Pa=zeros(U,L,length(Sam));
P=zeros(U,L);
uvar0 = 1;
uhat0 = 0;
%%
for i = 1 : length(Sam)    
  E = abs(Sam(1,i)*ones(U,L)-qhat).^2./qvar;
  Pa(:,:,i)=exp(-E);
  Pa = max(Pa,1e-200);
  P=P + Pa(:,:,i);
end                  
for i=1:4
  Pa(:,:,i)=Pa(:,:,i)./(P(:,:));
end  
%%
% Compute posterior mean and variance  
   for i=1 : length(Sam)
       xhat11 = Pa(:,:,i).*(Sam(1,i)*ones(U,L));
       xhat = xhat+xhat11;              
   end
   for i=1 : length(Sam)
       xvar11 = Pa(:,:,i).*(abs(Sam(1,i)*ones(U,L)-xhat).^2);
       xvar = xvar+xvar11; 
   end
   
   xvar_over_uvar0 = qvar./(uvar0+qvar);
    vall = (log(xvar_over_uvar0) + (1-xvar_over_uvar0) ...
        - abs(xhat-uhat0).^2./uvar0 );
end





